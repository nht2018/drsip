/*
 * Filename: Do not edit
 * Author: Hantao Nie (nht@pku.edu.cn)
 * Date: 2023-05-31 10:43:09
 * LastEditors: Hantao Nie (nht@pku.edu.cn)
 * LastEditTime: 2024-02-04 19:24:27
 * Description: 
 *      usage
 *      [X] = mexrotate_cone_r(Z, cone_size)
 *
 * Copyright (c) 2023, Hantao Nie, Peking University.
 */

// for rotated quadratic cone z \in K = {(x0, x_1, zbar) | 2 * x0 * x1 >= norm(zbar)}
// perform the transformation
//     (z0, z1, zbar) -> ( 1 / sqrt(2) *(z0 + z1), 1 / sqrt(2) * (z0 - z1), zbar)

#include "mex.h"
#include "omp.h"
#include <math.h>

#define EPS 1e-12 // numbers smaller than EPS are considered as zero
#define sqrt2 sqrt(2)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *Z, *X;
    double *cone_size_in;
    mwIndex *cone_size;
    mwIndex *cumsum_cone_size, *ind_head;
    mwIndex k; // number of blocks
    mwSize m_in, n_in, n;
    mwIndex i, j;

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

    // Get cone_size
    k = mxGetM(prhs[1]) * mxGetN(prhs[1]);
    cone_size_in = mxGetPr(prhs[1]);

    // transform cone_size to integer array
    cone_size = mxMalloc(k * sizeof(mwIndex));
#ifndef SERIAL
#pragma omp parallel for private(j) shared(cone_size, cone_size_in, k)
#endif
    for (j = 0; j < k; j++)
    {
        cone_size[j] = (mwIndex)cone_size_in[j];
    }

    // record the first index of each block
    ind_head = mxMalloc(k * sizeof(mwIndex));
    ind_head[0] = 0;
    for (j = 1; j < k; j++)
    {
        ind_head[j] = ind_head[j - 1] + cone_size[j - 1];
    }
    if (ind_head[k - 1] + cone_size[k - 1] != n)
        mexErrMsgTxt("cone_size does not match");

    // allocate output
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    X = mxGetDoubles(plhs[0]);

    // compute rotation
#ifndef SERIAL
#pragma omp parallel for private(j, i) shared(cone_size, ind_head, Z, X, k)
#endif
    for (j = 0; j < k; j++)
    {
        X[ind_head[j]] = (Z[ind_head[j]] + Z[ind_head[j] + 1]) / sqrt2;
        X[ind_head[j] + 1] = (Z[ind_head[j]] - Z[ind_head[j] + 1]) / sqrt2;
        for (i = 2; i < cone_size[j]; i++)
        {
            X[ind_head[j] + i] = Z[ind_head[j] + i];
        }
    }
    mxFree(cone_size);
    mxFree(ind_head);
}
