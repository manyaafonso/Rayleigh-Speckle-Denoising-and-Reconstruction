/*********************************************************************
 * Demo.cpp
 *
 * This file shows the basics of setting up a mex file to work with
 * Matlab.  This example shows how to use 2D matricies.  This may
 * 
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 *
 * For more information, see my site: www.shawnlankton.com
 * by: Shawn Lankton
 *
 ********************************************************************/
#include <matrix.h>
#include <mex.h>   

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void gradient3D(int dimx, int dimy, int dimz, double *img, double *dh, double *dv, double *dz);

void divergence3D(int dimx, int dimy, int dimz, double *x1, double *x2, double *x3, double *divx);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

//declare variables
    mxArray *a_in_m, *b_in_m, *c_in_m, *d_out_m;
    const mwSize *dims;
    double *a, *b, *c, *d;
    int dimx, dimy, dimz, numdims;
    int i,j,k;

//associate inputs
    a_in_m = mxDuplicateArray(prhs[0]);
    b_in_m = mxDuplicateArray(prhs[1]);
    c_in_m = mxDuplicateArray(prhs[2]);

//figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimz = (int)dims[0]; dimy = (int)dims[1]; dimx = (int)dims[2];


//associate outputs
    d_out_m = plhs[0] = mxCreateNumericArray( 3, dims, mxDOUBLE_CLASS, mxREAL );
    
//associate pointers
    a = mxGetPr(a_in_m);
    b = mxGetPr(b_in_m);
    c = mxGetPr(c_in_m);
    d = mxGetPr(d_out_m);
    
    divergence3D(dimx, dimy, dimz, a, b, c, d);

    return;
}




void divergence3D(int dimx, int dimy, int dimz, double *xr, double *xc, double *xs, double *divx) {
    
    int i,j,k;
    //do something
    for (i=0;i<dimx;i++) {
        
        for (j=0;j<dimy;j++) {
            
            for (k=0;k<dimz;k++) {
                
                if (k==0) { 
                    divx[i*dimy*dimz+j*dimz+k] = xr[i*dimy*dimz+j*dimz+k];
                }
                else {
                    if (k==dimz-1) {
                        divx[i*dimy*dimz+j*dimz+k] = -xr[i*dimy*dimz+j*dimz+dimz-2];
                    }
                    else {
                        divx[i*dimy*dimz+j*dimz+k] = xr[i*dimy*dimz+j*dimz+k]-xr[i*dimy*dimz+j*dimz+k-1];
                    }
                }
                
                if (j==0) { 
                    divx[i*dimy*dimz+j*dimz+k] = divx[i*dimy*dimz+j*dimz+k] + xc[i*dimy*dimz+j*dimz+k];
                }
                else {
                    if (j==dimy-1) {
                        divx[i*dimy*dimz+j*dimz+k] = divx[i*dimy*dimz+j*dimz+k] -xc[i*dimy*dimz+(dimy-2)*dimz+k];
                    }
                    else {
                        divx[i*dimy*dimz+j*dimz+k] = divx[i*dimy*dimz+j*dimz+k] + xc[i*dimy*dimz+j*dimz+k]-xc[i*dimy*dimz+(j-1)*dimz+k];
                    }
                }

                
                if (i==0) { 
                    divx[i*dimy*dimz+j*dimz+k] = divx[i*dimy*dimz+j*dimz+k] + xs[i*dimy*dimz+j*dimz+k];
                }
                else {
                    if (i==dimx-1) {
                        divx[i*dimy*dimz+j*dimz+k] = divx[i*dimy*dimz+j*dimz+k] -xs[(dimx-2)*dimy*dimz+j*dimz+k];
                    }
                    else {
                        divx[i*dimy*dimz+j*dimz+k] = divx[i*dimy*dimz+j*dimz+k] + xs[i*dimy*dimz+j*dimz+k]-xs[(i-1)*dimy*dimz+j*dimz+k];
                    }
                }

            }
        }
    }

}
