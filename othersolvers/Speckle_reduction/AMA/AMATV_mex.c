/*
% Paper   : Two-Level convex relaxed variational model for multiplicative denoising
% Authors : Myungjoo Kang, Sangwoon Yun, and Hyenkyun Woo(Corresponding Author)
% Journal : SIAM Journal of Imaging Sciences (submitted)
% Code written by Hyenkyun Woo (hyenkyun@gmail.com)
% Version : 2012-1-16-ver1.0
 */
        
#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include <math.h>
#include <time.h>

#define X(ix,iy) (ix)*iNy + (iy)
/*#define XX(ix,iy,i) (i)*((iNx+1)*(iNy+1)) + (ix)*(iNy+1) + (iy)*/
#define XX(ix,iy,i) (i)*((iNx+1)*(iNy+1)) + (iy)*(iNx+1) + (ix)
#define SQR(x) (x)*(x)
/*
float fast_log (float val)
{
   union {float f; int i;} t;
   t.f = val;   
   int* const exp_ptr = &t.i;
   int x              = *exp_ptr;
   const int log_2    = ((x >> 23) & 255) - 128;
   x &= ~(255 << 23);
   x += 127 << 23;
   *exp_ptr = x; 
   val = ((-1.0f/3) * t.f + 2) * t.f - 2.0f/3; 
    
   val = (val + log_2) * 0.69314718f;
   
   return (val);
}
*/
/****************************************/
extern void mexFunction(int iNbOut, mxArray *pmxOut[],
int iNbIn, const mxArray *pmxIn[])
{    
    float   *pfIm0, *pfu, *pfuOld, *pfuNew, *pfVecGeneralParameters, fmin, fmax,frho,fgrad,fhtemp1,fhtemp2;
    float   falpha, flambda, fGuij, fGuji, *pfz, *pfzNew,  *pfTV, fDivWD, fNormU, fTemp, fa;   
    float   *pfp, *pfpNew, fpij, fpji;
    int     iNy, iNx, iNdim, iDim[3],  iy, ix, i, j, iNyNx, iX, imodel;

    pfu    = mxGetData(pmxIn[0]);
    pfz    = mxGetData(pmxIn[1]);
    pfp    = mxGetData(pmxIn[2]);
    pfIm0  = mxGetData(pmxIn[3]);
  
    pfVecGeneralParameters = mxGetData(pmxIn[4]);
    
    iNy = (int) pfVecGeneralParameters[0];
    iNx = (int) pfVecGeneralParameters[1];
    falpha = pfVecGeneralParameters[2];       /*alpha*/
    flambda = pfVecGeneralParameters[3];      /*lambda*/
    fmin = pfVecGeneralParameters[4];
    fmax = pfVecGeneralParameters[5];
    imodel = (int) pfVecGeneralParameters[6];       /*model*/    
    
    iNdim = 2;
    iDim[0] = iNy;
    iDim[1] = iNx;
    
    pmxOut[0] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfuNew = mxGetData(pmxOut[0]);
        
    iNdim = 3;
    iDim[0] = iNx+1; /* extend to consider boundary */
    iDim[1] = iNy+1; 
    iDim[2] = 2;
    
    pmxOut[1] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfzNew = mxGetData(pmxOut[1]);
    
    pmxOut[2] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfpNew = mxGetData(pmxOut[2]);

    iNyNx = iNy*iNx;
    
        /* u */
       for(iy=0; iy< iNy; iy++)
       for(ix=0; ix< iNx; ix++)
       { 
            iX = X(ix,iy);
            fDivWD = 0.0;
           
            fpij = pfp[XX(ix+1,iy,  0)] - pfp[XX(ix,iy,0)];
            fpji = pfp[XX(ix,  iy+1,1)] - pfp[XX(ix,iy,1)];
            fDivWD = fpij + fpji +1;

            if(imodel==1 || imodel==2){ /* KL-DIV model */
                
               if(fDivWD == 0)
                   pfuNew[iX] = 1000*pfIm0[iX];
               else
                   pfuNew[iX] = pfIm0[iX]/fDivWD; 
               
            }else if(imodel==3 || imodel==4){ /* EXP-MIDAL model */
  
               if(fDivWD <= 0){
                   pfuNew[iX] = pfIm0[iX] + 10;   
               }else
                  /*pfuNew[iX] = pfIm0[iX] - fast_log(fDivWD);*/
                  pfuNew[iX] = pfIm0[iX] - log(fDivWD);
               
            }
            
            /*Projection of u*/
            if(pfuNew[iX] < fmin){
                pfuNew[iX] = fmin;
            }else if(pfuNew[iX] > fmax){
                pfuNew[iX] = fmax;
            }   
        
        }/*u : Gradient*/
      
        /* z */
        fa = 1/falpha;
        for(iy=0; iy< iNy+1; iy++)
        for(ix=0; ix< iNx+1; ix++)
        {
            fNormU = 0.0;
            
            fpij = pfp[XX(ix,iy,0)];
            fpji = pfp[XX(ix,iy,1)];
           
          if(ix>0 && ix<iNx && iy>0 && iy<iNy){ 
            fGuij = (pfuNew[X(ix,iy)]-pfuNew[X(ix-1,iy  )]);
            fGuji = (pfuNew[X(ix,iy)]-pfuNew[X(ix  ,iy-1)]);
          }else if(ix==0 && iy>0 && iy<iNy){
            fGuij = 0;
            fGuji = (pfuNew[X(ix,iy)]-pfuNew[X(ix  ,iy-1)]);  
          }else if(ix==iNx && iy>0 && iy<iNy){
            fGuij = 0;  
            fGuji = (pfuNew[X(ix-1,iy)]-pfuNew[X(ix-1,iy-1)]);
          }else if(iy==0 && ix>0 && ix<iNx){
            fGuij = (pfuNew[X(ix,iy)]-pfuNew[X(ix-1,iy)]);
            fGuji = 0;
          }else if(iy==iNy && ix>0 && ix<iNx){
            fGuij = (pfuNew[X(ix,iy-1)]-pfuNew[X(ix-1,iy-1)]);
            fGuji = 0;  
          }else if( (ix==0 && iy==0)||(ix==iNx && iy==0)||(ix==0 && iy==iNy)||(ix==iNx && iy==iNy) ){
            fGuij = 0;
            fGuji = 0;      
          }    
            
            pfzNew[XX(ix,iy,0)] = fGuij- fpij*fa;
            pfzNew[XX(ix,iy,1)] = fGuji- fpji*fa;
            
            fNormU = sqrt( SQR(pfzNew[XX(ix,iy,0)]) + SQR(pfzNew[XX(ix,iy,1)]) );
                                    
            
            if ( fNormU<flambda*fa ){
                
                pfzNew[XX(ix,iy,0)] = 0.0; /* z */
                pfzNew[XX(ix,iy,1)] = 0.0; 
                
                pfpNew[XX(ix,iy,0)] = pfp[XX(ix,iy,0)] - falpha*fGuij; /* p */
                pfpNew[XX(ix,iy,1)] = pfp[XX(ix,iy,1)] - falpha*fGuji;
                
            }else{
                
                fTemp = fNormU-flambda*fa; fTemp /= fNormU;
               
                pfzNew[XX(ix,iy,0)] *= fTemp; /* z  */
                pfzNew[XX(ix,iy,1)] *= fTemp; 
                    
                pfpNew[XX(ix,iy,0)] = pfp[XX(ix,iy,0)] - falpha*(fGuij - pfzNew[XX(ix,iy,0)]); /* p */
                pfpNew[XX(ix,iy,1)] = pfp[XX(ix,iy,1)] - falpha*(fGuji - pfzNew[XX(ix,iy,1)]); 
            } 
        }
}
