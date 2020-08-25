/*
 * permanent.c   Matlab MEX C-programming
 *   Compute the permanent of a matrix.  This routine seems to 
 * be efficient for sparse matrices.
 *
 * Copyright 2013-2016 Brian K. Butler
 */
/* Installation instructions:  First run "mex -setup" in MATLAB.
 * This routine compiles with either LCC, Microsoft Visual C++ 2008,
 * Express Edition, and MinGW64 Compiler (gcc).  Compile using 
 * "mex -v permanent.c" or "mex -O permanent.c" to optimize or even
 * "mex COPTIMFLAGS='-Ofast -Wall' permanent.c" for speed optimization
 * and warnings.  Now you are ready to run a test...
 *
 * Testing: Run the following in MATLAB
 *
    A=[1,0,0,0,0,0,0,0,1,0,0,1;...
        0,1,0,0,0,0,0,0,1,1,0,0;...
        0,0,1,0,0,0,0,0,0,1,1,0;...
        0,0,0,1,0,0,0,0,0,0,1,1;...
        0,0,0,0,1,0,0,0,1,1,1,0;...
        0,0,0,0,0,1,0,0,0,1,1,1;...
        0,0,0,0,0,0,1,0,1,0,1,1;...
        0,0,0,0,0,0,0,1,1,1,0,1;...
        0,0,0,0,1,1,0,0,1,0,0,0;...
        0,0,0,0,0,1,1,0,0,1,0,0;...
        0,0,0,0,0,0,1,1,0,0,1,0;...
        0,0,0,0,1,0,0,1,0,0,0,1];
    start=tic; p=permanent(A); t=toc(start)
 *
 * The resulting permanent should be 89.  t ~= 2.4e-04 on my machine.
 */
/* Changes from V1.0 to V1.1 Oct 28, 2015
 *   Added support for complex input matrices.
 *   Error checking for non-numeric and sparse format inputs.
 *   Return a permanent of 1 for 0x0 (empty matrix) input.
 * No changes to code from V1.1 to V1.2 
 * Changes from V1.2 to V1.3 Nov 14, 2016
 *   Added support for a matrix with more columns than rows.
 *   Switched the Laplace expansion from running along
 *   columns to running along rows.
 */
 
#include "matrix.h"
#include "mex.h"
#define A_IN (prhs[0])

static void   local_perm_cmplx(int currow, double *ppr, double *ppi);   
static double local_perm_real(int currow);
static void   prncrflgs( );
static int    m=0;      // Number of rows in matrix
static int    n=0;      // Number of cols in matrix
static double *a, *ai;  // pointer to matrix data
static int    *crflgs;  // col-removal flags

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *ppr, *ppi;
    
    if (nrhs != 1)
        mexErrMsgTxt("Function 'permanent' requires 1 input.\n");
    if (nlhs > 1)
        mexErrMsgTxt("Function 'permanent' only defined for 0 or 1 output(s).\n");
    if (!mxIsDouble(A_IN) || mxIsSparse(A_IN) || !mxIsNumeric(A_IN))
    {
        mexPrintf("Input is of type:\t%s\n",mxGetClassName(A_IN));
        mexErrMsgTxt("Function 'permanent' only uses double input (real or complex).\n");
    }
    m = mxGetM(A_IN);  // Number of rows in matrix
    n = mxGetN(A_IN);  // Number of cols in matrix

    //mexPrintf("num rows and cols %d %d\n", m, n);
	if (m > n)
    {
		mexErrMsgTxt("Matrix must be have # rows <= # columns.  Error inside permanent()\n");
    }

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxIsComplex(A_IN) ? mxCOMPLEX : mxREAL);
	ppr = mxGetPr(plhs[0]);

    if (m == 0)
	{
		*ppr = 1.0;  // 1 by definition.
	}
	else // m >= 1
	{
        crflgs = (int *)mxCalloc(n, sizeof(int));
        a = mxGetPr(A_IN);
      
        if (mxIsComplex(A_IN))
        {
            ai = mxGetPi(A_IN);
            ppi = mxGetPi(plhs[0]);
            local_perm_cmplx(0, ppr, ppi);            
        }
        else
            *ppr = local_perm_real(0);
        
        mxFree(crflgs);
	}
}

static double local_perm_real(int currow) {
// 'currow' indicates the current row, which we are to follow in the cofactor expansion.
    double p=0.0;
    int j;
    
    //mexPrintf("local_perm_real %d %d\t", m, currow);
    //mexPrintf("with flags\t "); prncrflgs( );

    if (m-currow == 1)  // down to a 1 row
    {
        for (j=0; j<n; j++) {
            if (crflgs[j])
                continue;
            
            p += a[currow+m*j];  // Sum not necessary for square matrix.  For non-square, need to sum remaining entries in final row.
        }
        //mexPrintf("1x1 p=%lf \n", p);
    }
    else {
        for (j=0; j<n; j++) {
            if (crflgs[j])
                continue;
            
            if (a[currow+m*j]) {
                crflgs[j]=1;
                p = p + a[currow+m*j] * local_perm_real(currow+1);
                crflgs[j]=0;
				//mexPrintf("weight %lf\n", a[currow+m*j]);
            }
        } // end for
    } //end else
    //mexPrintf("%dx%d perm=%lf \n", m-currow,n-currow,p);
    return (p);
}

static void local_perm_cmplx(int currow, double *pr, double *pi) {
// 'currow' indicates the current row, which we are to follow in the cofactor expansion.
    int j;
    double nr, ni;
    
    *pr=0.0;
    *pi=0.0;
    
    //mexPrintf("local_perm_cmplx %d %d\t", m, currow);
    //mexPrintf("with entries\t "); prncrflgs( );

    if (m-currow == 1)  // down to a 1x1
    {
        for (j=0; j<n; j++) {
            if (crflgs[j])
                continue;
            
            *pr += a[currow+m*j];
            *pi += ai[currow+m*j];
        }
        //mexPrintf("1x1 pr=%lf pi=%lf\n", *pr, *pi);
    }
    else {
        for (j=0; j<n; j++) {
            if (crflgs[j])
                continue;
            
            if (a[currow+m*j] || ai[currow+m*j]) {
                crflgs[j]=1;              
                local_perm_cmplx(currow+1, &nr, &ni);
                *pr = *pr + a[currow+m*j] * nr - ai[currow+m*j] * ni;
                *pi = *pi + ai[currow+m*j] * nr + a[currow+m*j] * ni;
                crflgs[j]=0;
            }
        } // end for
    } //end else
    //mexPrintf("%dx%d perm=%lf+j*%lf \n", m-currow,n-currow,*pr,*pi);
    return;
}

static void prncrflgs( ) 
{
    int j;
    
    for (j=0; j<n; j++) {
        if (crflgs[j])
            continue;       
        mexPrintf("%d  ", j);
    }
    mexPrintf("\n");
}
