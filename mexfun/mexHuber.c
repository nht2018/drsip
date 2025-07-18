/*
 * Filename: mexHuber.c
 * Author: Hantao Nie (nht@pku.edu.cn)
 * Date: 2023-05-31 10:43:09
 * LastEditors: Hantao Nie (nht@pku.edu.cn)
 * LastEditTime: 2024-03-15 10:01:28
 * Description: Huber function on R^n
 *      usage
 *      [X] = Huber(Z, epsilon)
 *
 * Copyright (c) 2023, Hantao Nie, Peking University.
 */


#define EPS 1e-12 // numbers smaller than EPS are considered as zero

#include "mex.h"
#include "omp.h"
#include <math.h>

inline double phi(double z, double epsilon)
{
    if (z <= 0)
        return 0;
    else if (z <= epsilon)
        return z * z / 2 / epsilon;
    else
        return z - epsilon / 2;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *Z, *X;
    double epsilon;
    mwSize m_in, n_in, n;

    if (nrhs != 2)
    {
        mexErrMsgTxt("Two inputs needed.");
        return;
    }

    // Get Z
    m_in = mxGetM(prhs[0]);
    n_in = mxGetN(prhs[0]);
    n = m_in * n_in;
    if (n_in != 1)
        mexErrMsgTxt("Z must be column vector.");

    Z = mxGetDoubles(prhs[0]);
    if (Z == NULL || mxIsSparse(prhs[0]))
    {
        mexErrMsgTxt("Z must be double");
        return;
    }


    // Get epsilon
    epsilon = mxGetScalar(prhs[1]);
    if (epsilon <= 0)
        mexErrMsgTxt("epsilon must be positive.");


    // allocate output
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    X = mxGetDoubles(plhs[0]);


    // compute 
    int j;
#ifndef SERIAL
#pragma omp parallel for private(j) shared(Z, X)
#endif
    for (j = 0; j < n; j++)
    {
        X[j] = phi(Z[j], epsilon);
    }
}
