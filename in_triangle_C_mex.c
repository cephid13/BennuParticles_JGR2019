///////////////////////////////////////////////////////////////////////////
// Test for point in triangle
// Author: Jay
///////////////////////////////////////////////////////////////////////////

#include    "mex.h"
#include	<math.h>
#include	<stdarg.h>
#include	<stdio.h>
#include	<stdlib.h>
#include    <string.h>

///////////////////
// -- Outputs -- //
///////////////////

double uv[2], *uv_Ptr;
int hit, *hit_Ptr;

//////////////////
// -- Inputs -- //
//////////////////

double *A, *B, *C, *P, *nhat;

/////////////////
// -- Locals -- //
/////////////////
double A_loc[3], B_loc[3], C_loc[3], P_loc[3], nhat_loc[3];
int mm;

/////////////////////
// -- Functions -- //
/////////////////////

void inverse ( double answer[], double x[], double y[], double n[] );
int in_triangle_3D( double *A, double *B, double *C, double *P, double *nhat, double *uv);

/***************************************************************
 * main program
 ***************************************************************/
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
        
{
    
    ///////////////////////////////////////////////////////////////////////
    
    /////////////////////////
    // -- Assign Inputs -- //
    /////////////////////////
    
    // - Vector and Matrix
    
    A       = mxGetPr(prhs[0]);       // 
    B       = mxGetPr(prhs[1]);       // 
    C       = mxGetPr(prhs[2]);       // 
    P       = mxGetPr(prhs[3]);       // 
    nhat    = mxGetPr(prhs[4]);       // 
    
    /*
    mexPrintf("facets[0] = %f\n",facets[0]);
    mexPrintf("facets[1] = %f\n",facets[1]);
    mexPrintf("facets[2] = %f\n",facets[2]);
    mexPrintf("facets[3] = %f\n",facets[3]);
    mexPrintf("facets[4] = %f\n",facets[4]);
    mexPrintf("facets[5] = %f\n",facets[5]);
    mexPrintf("verts[0] = %f\n",verts[0]);
    mexPrintf("verts[1] = %f\n",verts[1]);
    mexPrintf("verts[2] = %f\n",verts[2]);
    mexPrintf("verts[3] = %f\n",verts[3]);
     */
    
    ///////////////////////////////////////////////////////////////////////
    
    //////////////////////////
    // -- Assign Outputs -- //
    //////////////////////////
    
    plhs[0]  = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    hit_Ptr = (int *)mxGetData(plhs[0]);
    plhs[1]  = mxCreateDoubleMatrix(2, 1, mxREAL);
    uv_Ptr = mxGetPr(plhs[1]);
    
    ///////////////////////////////////////////////////////////////////////
    
    // copy data to local non const variables
    for (mm=0; mm<3; mm++) {
    	A_loc[mm] = A[mm];
    	B_loc[mm] = B[mm];
    	C_loc[mm] = C[mm];
    	P_loc[mm] = P[mm];
    	nhat_loc[mm] = nhat[mm];
    }
    
    // check if intersection is inside triangle
    //mexPrintf("P_loc[0] = %f\n",P_loc[0]);
    hit =  in_triangle_3D( A_loc, B_loc, C_loc, P_loc, nhat_loc, uv);
    //mexPrintf("hit = %d\n",hit);
    //mexPrintf("uv[0] = %f\n",uv[0]);
    //mexPrintf("uv[1] = %f\n",uv[1]);
    
    // assign outputs
    hit_Ptr[0] = hit;
    uv_Ptr[0] = uv[0];
    uv_Ptr[1] = uv[1];
    
    return;
    
} // For main

///////////////////////////////////////////////////////////////////////////

////// SUBFUNCTIONS //////

void inverse ( double answer[6], double x[3], double y[3], double n[3] )
{
    double denom;
    
    denom = x[0]*(y[1]*n[2]-n[1]*y[2]) - y[0]*(x[1]*n[2]-n[1]*x[2]) + n[0]*(x[1]*y[2]-y[1]*x[2]);
    
    answer[0] = (y[1]*n[2]-n[1]*y[2])/denom;
    answer[1] = (y[2]*n[0]-n[2]*y[0])/denom;
    answer[2] = (y[0]*n[1]-n[0]*y[1])/denom;
    answer[3] = (x[2]*n[1]-n[2]*x[1])/denom;
    answer[4] = (x[0]*n[2]-n[0]*x[2])/denom;
    answer[5] = (x[1]*n[0]-n[1]*x[0])/denom;
    //    answer[6] = (y[2]*x[1]-x[2]*y[1])/denom;
    //    answer[7] = (y[0]*x[2]-x[0]*y[2])/denom;
    //    answer[8] = (y[1]*x[0]-x[1]*y[0])/denom;
    
    return;
}

int in_triangle_3D( double *A, double *B, double *C, double *P, double *nhat, double *uv)
{
    // In triangle stuff
    double v0[3], v1[3], v2[3], answer[6];
    int i;
    
    for (i=0;  i<3; i++){
        v0[i] = C[i] - A[i];
        v1[i] = B[i] - A[i];
        v2[i] = P[i] - A[i];
    }
    
    //mexPrintf("v2[0] = %f\n",v2[0]);
    
    inverse( answer, v1, v0, nhat );
    
    uv[0] = answer[0]*v2[0] + answer[1]*v2[1] + answer[2]*v2[2];
    uv[1] = answer[3]*v2[0] + answer[4]*v2[1] + answer[5]*v2[2];
    
    return ((uv[0] >= 0) && (uv[1] >= 0) && (uv[0] + uv[1] <= 1));
}


///////////////////////////////////////////////////////////////////////////