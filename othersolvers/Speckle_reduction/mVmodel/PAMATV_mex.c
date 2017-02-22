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

float fast_pow(float base,int exp)
{
    float result = base;
    
    switch(exp){
    case 3:  base *= base; break;   
    case 4:  base *= base;    base *= result;  break;    
    case 5:  base *= base;    base *= base;    break;
    case 9:  base *= base;    base *= base;  base *= base;  break;
    case 17: base *= base;    base *= base;  base *= base; base *= base;  break;
    case 33: base *= base;    base *= base;  base *= base; base *= base; base *= base;  break;
    case 65: base *= base;    base *= base;  base *= base; base *= base; base *= base; base *= base; break;   
    default: 
        result = pow(base,exp);  
    break;
    }
    
    result *= base;
    
    return result;
}



/****************************************/
extern void mexFunction(int iNbOut, mxArray *pmxOut[],
int iNbIn, const mxArray *pmxIn[])
{    
    float   *pfIm0, *pfu, *pfuOld, *pfuNew, *pfVecGeneralParameters, fmin, fmax,frho,fgrad,fhtemp1,fhtemp2;
    float   falpha, flambda, fQQ, fGuij, fGuji, *pfz, *pfzNew,  *pfTV, fDivWD, fNormU, fTemp;   
    float   *pfp, *pfpNew, fpij, fpji, fa, fla;
    int     iNy, iNx, iNdim, iDim[3], iDual, iy, ix, i, j, iNyNx, iX, igsq, igsq1;

    pfu    = mxGetData(pmxIn[0]);
    pfz    = mxGetData(pmxIn[1]);
    pfp    = mxGetData(pmxIn[2]);
    pfIm0  = mxGetData(pmxIn[3]);
  
    pfVecGeneralParameters = mxGetData(pmxIn[4]);
    
    iNy = (int) pfVecGeneralParameters[0];
    iNx = (int) pfVecGeneralParameters[1];
    falpha = pfVecGeneralParameters[2];       /*alpha*/
    flambda = pfVecGeneralParameters[3];      /*lambda*/
    igsq = (int) pfVecGeneralParameters[4];
    fmin = pfVecGeneralParameters[5];
    fmax = pfVecGeneralParameters[6];
    frho = pfVecGeneralParameters[7];          /*rho*/    

    iNdim = 2;
    iDim[0] = iNy;
    iDim[1] = iNx;
    
    pmxOut[0] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfuNew = mxGetData(pmxOut[0]);
        
    iNdim = 3;
    iDim[0] = iNx+1; /*iNy+1*/ /* extend to consider boundary */
    iDim[1] = iNy+1; /*iNx+1*/
    iDim[2] = 2;
    
    pmxOut[1] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfzNew = mxGetData(pmxOut[1]);
    
    pmxOut[2] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfpNew = mxGetData(pmxOut[2]);
         
    iNyNx = iNy*iNx;
    
        /* u */
       igsq1 = igsq+1;
       for(iy=0; iy< iNy; iy++)
       for(ix=0; ix< iNx; ix++)
       { 
            iX = X(ix,iy);
            fDivWD = 0.0;
           
            fpij = pfp[XX(ix+1,iy,  0)] - pfp[XX(ix,iy,0)];
            fpji = pfp[XX(ix,  iy+1,1)] - pfp[XX(ix,iy,1)];
            fDivWD = fpij + fpji;

            fgrad = (igsq/pfu[iX])- pfIm0[iX]/fast_pow(pfu[iX],igsq1);            
            pfuNew[iX] = pfu[iX] - frho*(fgrad + fDivWD);  
            
            /*Projection of u*/
            if(pfuNew[iX] < fmin){
                pfuNew[iX] = fmin;
            }else if(pfuNew[iX] > fmax){
                pfuNew[iX] = fmax;
            }   
        
        }/*u : Gradient*/
      
        /* z */
        fa = 1/falpha;
        fla = flambda*fa;
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
                                    
            
            if ( fNormU<fla ){
                
                pfzNew[XX(ix,iy,0)] = 0.0; /* z */
                pfzNew[XX(ix,iy,1)] = 0.0; 
                
                pfpNew[XX(ix,iy,0)] = pfp[XX(ix,iy,0)] - falpha*fGuij; /* p */
                pfpNew[XX(ix,iy,1)] = pfp[XX(ix,iy,1)] - falpha*fGuji;
                
            }else{
                
                fTemp = fNormU-fla; fTemp /= fNormU;
               
                pfzNew[XX(ix,iy,0)] *= fTemp; /* z  */
                pfzNew[XX(ix,iy,1)] *= fTemp; 
                    
                pfpNew[XX(ix,iy,0)] = pfp[XX(ix,iy,0)] - falpha*(fGuij - pfzNew[XX(ix,iy,0)]); /* p */
                pfpNew[XX(ix,iy,1)] = pfp[XX(ix,iy,1)] - falpha*(fGuji - pfzNew[XX(ix,iy,1)]); 
            } 
        }
}
