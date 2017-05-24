/* 	
 function m_regul_totvar_2d

  [fx, gx] = m_regul_totvar_2d(x, w1, w2, eps, flags);
                                      
  MexFile to call from Matlab the rgl_tv2d function of the Totvar library.
 
  Compilation :
        mex m_regul_totvar_2d.c totvar.o
  
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
    double w1, w2, eps;
    mwSize n1, n2, nDims; 
    const mwSize *pDims;
    const mwSize pDims2[] = {1,1};
	unsigned int flags; 
	double *x;

    double *fx;
    double *gx; 

    /* Control of the number of inputs and outputs */ 
    if (nrhs > 5)
        mexErrMsgTxt("5 input argument required."); 
    else if (nlhs1 !=2) 
        mexErrMsgTxt("2 output argument required.");
    
    /* Control of the inputs */
	if (!mxIsNumeric(prhs[0]) || !mxIsDouble(prhs[0])) 
        mexErrMsgTxt("The first input must be numerical double array.");
    if (!mxIsNumeric(prhs[1]) || !mxIsScalar(prhs[1])
        || !mxIsNumeric(prhs[2]) || !mxIsScalar(prhs[2]) 
		|| !mxIsNumeric(prhs[3]) || !mxIsScalar(prhs[3])
        || !mxIsNumeric(prhs[4]) || !mxIsScalar(prhs[4]) ) 
        mexErrMsgTxt("The 4 last input must be scalar."); 
    
    if ( (double )(unsigned int )mxGetScalar(prhs[4]) != mxGetScalar(prhs[4]) )
        mexErrMsgTxt("The last input flags must be an unsigned integer.");
    
    /* Get the inputs */
    w1        = (double )mxGetScalar(prhs[1]);
    w2        = (double )mxGetScalar(prhs[2]);
    eps       = (double )mxGetScalar(prhs[3]);
    flags     = (unsigned int )mxGetScalar(prhs[4]);
    flags    |= RGL_INTEGRATE_GRADIENT;
		
    nDims = mxGetNumberOfDimensions(prhs[0]);
    pDims = mxGetDimensions(prhs[0]);
    x = mxGetPr(prhs[0]);	
    
    n1        = pDims[0];
    n2        = pDims[1];

    /* Display the inputs */
    /* mexPrintf("n1 = %d, n2 = %d, w1 = %f, w2 = %f, eps = %f, flags = %d\n", n1, n2, w1, w2, eps, flags); */
    
    /* Create the output arguments */
    /* Argument fx */
    plhs[0] = mxCreateNumericArray(2, pDims2, mxDOUBLE_CLASS, mxREAL);
    fx = mxGetPr(plhs[0]);
    
    /* Argument gx */
    plhs[1] = mxCreateNumericArray(nDims, pDims, mxDOUBLE_CLASS, mxREAL);
    gx = mxGetPr(plhs[1]);
   
    /* Call to the rgl_tv3d function */
    *fx = rgl_tv2d(x,n1,n2,w1,w2,eps,gx,flags); 
}
