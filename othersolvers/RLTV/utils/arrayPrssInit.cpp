/* 
 * =============================================================
 * arrayPrssInit.cpp
 * Initialize denoising algorithm (implemented in C++)
 * =============================================================
 */
#include <math.h>
#include "mex.h"

int c3(int &l, int &c, int &p, int dim[3]) //access given matrix coordinates (x,y,z) in C
{
    int nl=dim[0], nc=dim[1], np=dim[2], npc=nl*nc;
    return(npc*p+c*nl+l);
}
double Nv(int &i, int &j, int &k, int dim[3], double *g)
{
    double res;
    int i2, j2, k2;
    i2=(i<dim[0]-1?i+1:dim[0]-2); j2=(j<dim[1]-1?j+1:dim[1]-2); k2=(k<dim[2]-1?k+1:dim[2]-1-(dim[2]>1));
    res = 3/g[c3(i,j,k,dim)] + 1/g[c3(i2,j,k,dim)] + 1/g[c3(i,j2,k,dim)] + 1/g[c3(i,j,k2,dim)];
    return(res);
}
double meanX(int &i, int &j, int &k, int dim[3], double *g, double *x)
{
    int i1, j1, k1, i2, j2, k2;
    
    i1=(i>0?i-1:1);   j1=(j>0?j-1:1);   k1=(k>0?k-1:dim[2]>1);
    i2=(i<dim[0]-1?i+1:dim[0]-2); j2=(j<dim[1]-1?j+1:dim[1]-2); k2=(k<dim[2]-1?k+1:dim[2]-1-(dim[2]>1)); 
    return(1/Nv(i,j,k,dim,g) * ((1/g[c3(i,j,k,dim)])*(x[c3(i1,j,k,dim)]+x[c3(i,j1,k,dim)]+x[c3(i,j,k1,dim)]) + (1/g[c3(i2,j,k,dim)])*x[c3(i2,j,k,dim)] + (1/g[c3(i,j2,k,dim)])*x[c3(i,j2,k,dim)] + (1/g[c3(i,j,k2,dim)])*x[c3(i,j,k2,dim)]));
}
void gradient(double *x, double *g, int dim[3])
{
    double dx, dy, dz;
    int i1, j1, k1; 
  

    for (int k = 0; k < dim[2]; k++)            //sweeping planes
        for (int i = 0; i < dim[0]; i++)        // sweeping lines
            for (int j = 0; j < dim[1]; j++)    //sweeping column
            {
                i1=(i>0?i-1:1);   j1=(j>0?j-1:1);   k1=(k>0?k-1:dim[2]>1);
                
                dx=x[c3(i,j,k,dim)]-x[c3(i,j1,k,dim)];
                dy=x[c3(i,j,k,dim)]-x[c3(i1,j,k,dim)];
                dz=x[c3(i,j,k,dim)]-x[c3(i,j,k1,dim)];
            
                g[c3(i,j,k,dim)]= sqrt(dx*dx + dy*dy + dz*dz + 0.001);              
            }
}
void arrayPrssInit(double *x, double *y, double *x_n, int dim[3], double alpha)
{
   double *g        = new double[dim[0]*dim[1]*dim[2]];
   double E, nablaE, nabla2E;
   
   gradient(x,g,dim);
      
   int i1, j1, k1;
   double dx2, dy2, dz2;
   
   
   for (int k = 0; k < dim[2]; k++)            // sweeping planes
       for (int i = 0; i < dim[0]; i++)        // sweeping lines
           for (int j = 0; j < dim[1]; j++)    // sweeping column
           {
               i1=(i>0?i-1:1);   j1=(j>0?j-1:1);   k1=(k>0?k-1:(dim[2]>1));               

               dx2=pow(x[c3(i,j,k,dim)]-x[c3(i,j1,k,dim)],2);
               dy2=pow(x[c3(i,j,k,dim)]-x[c3(i1,j,k,dim)],2);
               dz2=pow(x[c3(i,j,k,dim)]-x[c3(i,j,k1,dim)],2);
               
            
               // Compute E, nablaE and nabla2E
               E = 0.5*pow(y[c3(i,j,k,dim)],2)*exp(-x[c3(i,j,k,dim)]) + x[c3(i,j,k,dim)] + 
               alpha*(1/g[c3(i,j,k,dim)])*(dx2+dy2+dz2); 
               
               nablaE = -0.5*pow(y[c3(i,j,k,dim)],2)*exp(-x[c3(i,j,k,dim)]) + 1 + 
               2*alpha*Nv(i,j,k,dim,g)*(x[c3(i,j,k,dim)] - meanX(i,j,k,dim,g,x));         
               
               nabla2E = 0.5*pow(y[c3(i,j,k,dim)],2)*exp(-x[c3(i,j,k,dim)]) + 2*alpha*Nv(i,j,k,dim,g);
               
               x_n[c3(i,j,k,dim)] = x[c3(i,j,k,dim)] - 0.2*(nablaE/nabla2E);  
            }
    delete [] g;
}

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
//nlhs    The number of expected output mxArrays
//plhs    Array of pointers to the expected output mxArrays
//nrhs    The number of input mxArrays
//prhs    Array of pointers to the input mxArrays. 
//
//
// Get the number of columns                            N = mxGetN(prhs[0]) 
// Get the number of lines                              M = mxGetM(prhs[0]) 
// Get the scalar input x.                              x = mxGetScalar(prhs[0])
// Create a pointer to the input matrix y.              y = mxGetPr(prhs[1]);
// Set the output pointer to the output matrix.         plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
// Create a C pointer to a copy of the output matrix.   z = mxGetPr(plhs[0]);
// Check to make sure the first input argument is a scalar.
// if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetN(prhs[0])*mxGetM(prhs[0]) != 1)
//      mexErrMsgTxt("Input x must be a scalar.");
// ----------------------------------------------------    
    
  double *x, *y, alpha;             //Inputs
  double *x_n;                      //Outputs
  
  
  /*  Check for proper number of arguments. */
  if (nrhs != 3) mexErrMsgTxt("Three inputs required.");
  if (nlhs != 1) mexErrMsgTxt("One output required.");
  
  
  /* INPUTS */
  
  x     = mxGetPr(prhs[0]); // Create a pointer to the input matrix x.
  y     = mxGetPr(prhs[1]); // Create a pointer to the input matrix y.

  // Get number and dimensions of x
  int nDim=0, dim[3]={1,1,1};
  const int *n;
  nDim  = mxGetNumberOfDimensions(prhs[0]);
  n     = mxGetDimensions(prhs[0]);
  alpha = mxGetScalar(prhs[2]);
  
  for(int i=0;i<nDim;i++) dim[i]=n[i];
  
  plhs[0] = mxCreateDoubleMatrix(dim[0],dim[1]*dim[2], mxREAL); //Set the output pointer to the output matrix.  
  x_n      = mxGetPr(plhs[0]); //Create a C pointer to a copy of the output matrix.
   
  arrayPrssInit(x,y,x_n,dim,alpha); // Call the Cpp subroutine.
}