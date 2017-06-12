/* 	
 function m_regul_totvar_3d_complex

  [fx, gxr, gxi] = m_regul_totvar_3d_complex(xr, xi, w1, w2, w3, eps, flags);
                                      
  MexFile to call from Matlab the rgl_tv3d function of the Totvar library
  which implements Definitions for relaxed Total Variation (TV)
  for images with complex values.
 
  Compilation :
        mex m_regul_totvar_3d_complex.c ../totvar.o (or totvar.obj)
  
  Juillet 2015
  CEA-LETI
  Fabien Momey
*/

#include "matrix.h"
#include "mex.h"
#include "../totvar.h"

void mexFunction( int nlhs1, mxArray *plhs[],    
                  int nrhs, const mxArray *prhs[] )
{ 
    double w1, w2, w3, eps;
    mwSize n1, n2, n3, nDims; 
    const mwSize *pDims;
    const mwSize pDims2[] = {1,1};
	unsigned int flags; 
	double *xr;
    double *xi;

    double *fx;
    double *gxr; 
    double *gxi; 

    /* Control of the number of inputs and outputs */ 
    if (nrhs > 7)
        mexErrMsgTxt("7 input argument required."); 
    else if (nlhs1 !=3) 
        mexErrMsgTxt("3 output argument required.");
    
    /* Control of the inputs */
	if (!mxIsNumeric(prhs[0]) || !mxIsDouble(prhs[0])) 
        mexErrMsgTxt("The first input must be numerical double array.");
    if (!mxIsNumeric(prhs[1]) || !mxIsDouble(prhs[1])) 
        mexErrMsgTxt("The second input must be numerical double array.");
    if (!mxIsNumeric(prhs[2]) || !mxIsScalar(prhs[2])
        || !mxIsNumeric(prhs[3]) || !mxIsScalar(prhs[3]) 
		|| !mxIsNumeric(prhs[4]) || !mxIsScalar(prhs[4])
        || !mxIsNumeric(prhs[5]) || !mxIsScalar(prhs[5]) 
		|| !mxIsNumeric(prhs[6]) || !mxIsScalar(prhs[6]) ) 
        mexErrMsgTxt("The 5 last input must be scalar."); 
    
    if ( (double )(unsigned int )mxGetScalar(prhs[6]) != mxGetScalar(prhs[6]) )
        mexErrMsgTxt("The last input flags must be an unsigned integer.");
    
    /* Get the inputs */
    w1        = (double )mxGetScalar(prhs[2]);
    w2        = (double )mxGetScalar(prhs[3]);
	w3        = (double )mxGetScalar(prhs[4]);
    eps       = (double )mxGetScalar(prhs[5]);
    flags     = (unsigned int )mxGetScalar(prhs[6]);
    flags    |= RGL_INTEGRATE_GRADIENT;
		
	nDims = mxGetNumberOfDimensions(prhs[0]);
    pDims = mxGetDimensions(prhs[0]);
    xr = mxGetPr(prhs[0]);
    
    xi = mxGetPr(prhs[1]);
    
    n1        = pDims[0];
    n2        = pDims[1];
	n3        = pDims[2];

    /* Display the inputs */
    /* mexPrintf("n1 = %d, n2 = %d, n3 = %d, w1 = %f, w2 = %f, w3 = %f, eps = %f, flags = %d\n", n1, n2, n3, w1, w2, w3, eps, flags); */
    
    /* Create the output arguments */
    /* Argument fx */
    plhs[0] = mxCreateNumericArray(2, pDims2, mxDOUBLE_CLASS, mxREAL);
    fx = mxGetPr(plhs[0]);
    
    /* Argument gxr */
    plhs[1] = mxCreateNumericArray(nDims, pDims, mxDOUBLE_CLASS, mxREAL);
    gxr = mxGetPr(plhs[1]);
    
    /* Argument gxi */
    plhs[2] = mxCreateNumericArray(nDims, pDims, mxDOUBLE_CLASS, mxREAL);
    gxi = mxGetPr(plhs[2]);
   
    /* Call to the rgl_tv3d function */
    *fx = rgl_tv3d_complex(xr,xi,n1,n2,n3,w1,w2,w3,eps,gxr,gxi,flags); 
}
