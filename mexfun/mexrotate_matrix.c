/*
 * Filename: Do not edit
 * Author: Hantao Nie (nht@pku.edu.cn)
 * Date: 2023-05-31 10:43:09
 * LastEditors: Hantao Nie (nht@pku.edu.cn)
 * LastEditTime: 2023-10-27 13:39:47
 * Description: the proximal onto the second order cone
 *      usage
 *      [T] = mexrotate_matrix(Z, cone_size)
 *
 * Copyright (c) 2023, Hantao Nie, Peking University.
 */

// for rotated quadratic cone z \in K = {(x0, x_1, zbar) | 2 * x0 * x1 >= norm(zbar)}
// obtain the rotation matrix
// T = [1 / sqrt(2), 1 / sqrt(2) , 0;
//      1 / sqrt(2), -1 / sqrt(2), 0;
//       0, 0, I]

#include "mex.h"
#include "omp.h"
#include <math.h>

#define EPS 1e-12 // numbers smaller than EPS are considered as zero
#define sqrt2 sqrt(2)
#define sqrt05 1 / sqrt2

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *cone_size_in;
    mwIndex *cone_size;
    mwIndex *cumsum_cone_size, *ind_head;
    mwIndex k; // number of blocks
    mwSize m_in, n_in, n;
    mwIndex i, j;
    mwIndex *ir, *jc;
    double *pr;

    if (nrhs != 1)
    {
        mexErrMsgTxt("One input needed.");
        return;
    }

    // Get cone_size
    k = mxGetM(prhs[0]) * mxGetN(prhs[0]);
    cone_size_in = mxGetPr(prhs[0]);

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
    n = ind_head[k - 1] + cone_size[k - 1];

    // allocate output
    plhs[0] = mxCreateSparse(n, n, n + 2 * k, mxREAL);
    ir = mxGetIr(plhs[0]);
    jc = mxGetJc(plhs[0]);
    pr = mxGetPr(plhs[0]);

    int index;
#ifndef SERIAL
#pragma omp parallel for private(j, i, index) shared(ir, jc, pr, k, ind_head) schedule(dynamic, 1)
#endif
    for (j = 0; j < k; j++)
    {
        index = ind_head[j] + 2 * j;
        ir[index] = ind_head[j];
        ir[index + 1] = ind_head[j] + 1;
        ir[index + 2] = ind_head[j];
        ir[index + 3] = ind_head[j] + 1;
        pr[index] = sqrt05;
        pr[index + 1] = sqrt05;
        pr[index + 2] = sqrt05;
        pr[index + 3] = -sqrt05;
        jc[ind_head[j]] = index;
        jc[ind_head[j] + 1] = index + 2;
        index = index + 2;
        for (i = 2; i < cone_size[j]; i++)
        {
            ir[index + i] = ind_head[j] + i;
            pr[index + i] = 1;
            jc[ind_head[j] + i] = index + i;
        }
    }
    jc[n] = n + 2 * k;

    mxFree(cone_size);
    mxFree(ind_head);
}
