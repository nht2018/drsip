/*
 * Filename: Do not edit
 * Author: Hantao Nie (nht@pku.edu.cn)
 * Date: 2023-05-31 10:43:09
 * LastEditors: Hantao Nie (nht@pku.edu.cn)
 * LastEditTime: 2023-10-27 13:39:59
 * Description: the projection onto the second order cone
 *      usage
 *      [X] = mexprojection_cone_q(Z, cone_size)
 *      [X, eigval, Omega1, Omega2, U1, U2, shift] = mexprojection_cone_q(Z, cone_size)
 *
 * Copyright (c) 2023, Hantao Nie, Peking University.
 */

// compute the projection onto the second order cone also its Jacobian
// For a single quadratic cone K = {(x0, xbar) | x0 >= norm(xbar)}
// the spectral decomposition is given by
// z = lambda1 * c1 * c1' + lambda2 * c2 * c2'
// where lambda1 = z0 + norm(zbar), lambda2 = z0 + norm(zbar),
// c1 = 0.5 * [1; zbar / norm(zbar)], c2 = 0.5 * [1; -zbar / norm(zbar)]
// the proximal mapping is given by
//      prox(x) = max(lambda1, 0) * c1 + max(lambda2, 0) * c2
// the derivative of prox(x) is given by
//      Dprox(x) = d * I + U * Omega * U'
// with d = (max(lambda1, 0) - max(lambda2, 0)) / (2 * norm(zbar))
//      U = sqrt(2) * [c1, c2]
//      Omega = diag(Dmax(lambda1, 0) - d, Dmax(lambda2, 0) - d )

// Assume the input Z has size [n, 1] and consists of k subblocks. The output of each variable is as follows:
// X: projection of Z onto the cone, i.e., X = prox(Z) , size = [n, 1]
// eigval: [lambda1, lambda2], size = [k, 2]
// Omega1: size = [k, 1]
// Omega2: size = [k, 1]
// U1: sqrt(2) * c1, size = [n, 1]
// U2: sqrt(2) * c2, size = [n, 1]
// shift: d, size = [k, 1]

#define EPS 1e-12 // numbers smaller than EPS are considered as zero

#include "mex.h"
#include "omp.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *Z, *X;
    double *cone_size_in;
    mwIndex *cone_size;
    mwIndex *cumsum_cone_size, *ind_head;
    mwIndex k; // number of blocks
    mwSize m_in, n_in, n;
    mwIndex i, j;
    double *Omega1, *Omega2, *eigval, *U1, *U2, *shift;

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
        mexErrMsgTxt("cone_size does not match X");

    // allocate output
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    X = mxGetDoubles(plhs[0]);
    if (nlhs > 1)
    {
        plhs[1] = mxCreateDoubleMatrix(k, 2, mxREAL);
        eigval = mxGetDoubles(plhs[1]);
        plhs[2] = mxCreateDoubleMatrix(k, 1, mxREAL);
        Omega1 = mxGetDoubles(plhs[2]);
        plhs[3] = mxCreateDoubleMatrix(k, 1, mxREAL);
        Omega2 = mxGetDoubles(plhs[3]);
        plhs[4] = mxCreateDoubleMatrix(n, 1, mxREAL);
        U1 = mxGetDoubles(plhs[4]);
        plhs[5] = mxCreateDoubleMatrix(n, 1, mxREAL);
        U2 = mxGetDoubles(plhs[5]);
        plhs[6] = mxCreateDoubleMatrix(k, 1, mxREAL);
        shift = mxGetDoubles(plhs[6]);
    }

    // compute projection
    double z0, detz, norm2zbar, normzbar;
#ifndef SERIAL
#pragma omp parallel for private(j, i, z0, detz, norm2zbar, normzbar) shared(cone_size, ind_head, k, Z, X, eigval, Omega1, Omega2, U1, U2, shift)
#endif
    for (j = 0; j < k; j++)
    {
        z0 = Z[ind_head[j]];
        norm2zbar = 0;
        for (i = 1; i < cone_size[j]; i++)
        {
            norm2zbar += Z[ind_head[j] + i] * Z[ind_head[j] + i];
        }
        detz = z0 * z0 - norm2zbar;
        normzbar = sqrt(norm2zbar);

        if (detz >= 0 && z0 > 0)
        { // in the interior of the cone
            for (i = 0; i < cone_size[j]; i++)
                X[ind_head[j] + i] = Z[ind_head[j] + i];
        }
        else if (detz >= 0 && z0 <= 0)
        { // in the interior of the mirror of the cone
            for (i = 0; i < cone_size[j]; i++)
                X[ind_head[j] + i] = 0;
        }
        else
        { // otherwise
            X[ind_head[j]] = 0.5 * (z0 + normzbar);
            for (i = 1; i < cone_size[j]; i++)
                X[ind_head[j] + i] = 0.5 * (z0 / normzbar + 1) * Z[ind_head[j] + i];
        }

        // compute Jacobian of the projection
        if (nlhs > 1)
        {
            eigval[j] = z0 + normzbar;
            eigval[k + j] = z0 - normzbar; // eigval[k+j] in c is eigval(j, 2) in matlab
            U1[ind_head[j]] = sqrt(0.5);
            U2[ind_head[j]] = sqrt(0.5);
            for (i = 1; i < cone_size[j]; i++)
            {
                if (fabs(normzbar) < EPS)
                    U1[ind_head[j] + i] = 0.;
                else
                    U1[ind_head[j] + i] = sqrt(0.5) * Z[ind_head[j] + i] / normzbar;
                U2[ind_head[j] + i] = -U1[ind_head[j] + i];
            }

            if (detz >= 0 && z0 > 0)
            { // in the interior of the cone
                shift[j] = 1;
                Omega1[j] = 0;
                Omega2[j] = 0;
            }
            else if (detz >= 0 && z0 <= 0)
            { // in the interior of the mirror of the cone
                shift[j] = 0;
                Omega1[j] = 0;
                Omega2[j] = 0;
            }
            else
            { // otherwise
                shift[j] = 0.5 * (1 + z0 / normzbar);
                Omega1[j] = 1 - shift[j];
                Omega2[j] = -shift[j];
            }
        }
    }
    mxFree(cone_size);
    mxFree(ind_head);
}
