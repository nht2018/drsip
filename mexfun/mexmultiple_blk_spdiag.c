
/*

input arguments
    A: a sparse matrix with size [m, n]
    X: a double vector with size [n, 1] or [1, n]
    block_size: a integer vector with length k, i.e., block_size = [q_1 ,..., q_k]
        where q_1 + ... + q_k == n

output
    A * blk_spdiag(X, block_size)


compile instruction
    this code need openmp
    on Linux
        mex -c CFLAGS='-fopenmp' blk_spdiag.c
        mex -L/path/to/matlab/sys/os/glnxa64 -liomp5 -lpthread blk_spdiag.o
    on macOS,
        curl -O https://mac.r-project.org/openmp/openmp-13.0.0-darwin21-Release.tar.gz
        sudo tar fvxz openmp-13.0.0-darwin21-Release.tar.gz -C /
        mex -R2018a CXX_FLAGS="-Xclang -fopenmp" LDFLAGS="$LDFLAGS -lomp" CXXOPTIMFLAGS="$CXXOPTIMFLAGS -Xclang -fopenmp" -I/usr/local/include blk_spdiag.c
        see https://ww2.mathworks.cn/matlabcentral/answers/1761950-m1-mac-compile-mex-file-with-openmp
        But it seems not work on my laptop; The executable file is still single-thread.

*/

#include "mex.h"
#include "omp.h"
#include "string.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *X;
    double *cone_size_in;
    mwIndex *block_size;
    mwIndex *cumsum_cone_size;
    mwIndex k; // number of blocks
    mwSize m_in, n_in, n;
    mwIndex i, j;

    mwIndex *A_ir, *A_rc;
    double *A_pr;

    mwIndex *ir, *jc;
    double *pr;

    // Check number of input arguments
    if (nrhs != 3)
        mexErrMsgTxt("Three  input arguments required.");

    // Get A
    A_ir = mxGetIr(prhs[0]);
    A_rc = mxGetJc(prhs[0]);
    A_pr = mxGetPr(prhs[0]);


    // Get X
    m_in = mxGetM(prhs[1]);
    n_in = mxGetN(prhs[1]);
    n = m_in * n_in;
    if (!(m_in == 1 || n_in == 1))
        mexErrMsgTxt("X must be vector.");

    X = mxGetPr(prhs[0]);

    // Get block_size
    k = mxGetM(prhs[2]) * mxGetN(prhs[2]);
    cone_size_in = mxGetPr(prhs[2]);

    block_size = mxMalloc(k * sizeof(mwIndex));
#ifndef SERIAL
#pragma omp parallel for private(j) shared(block_size, cone_size_in, k)
#endif
    for (j = 0; j < k; j++)
    {
        block_size[j] = (mwIndex)cone_size_in[j];
    }
    cumsum_cone_size = mxMalloc((k + 1) * sizeof(mwIndex));
    cumsum_cone_size[0] = 0;
    for (j = 0; j < k; j++)
    {
        cumsum_cone_size[j + 1] = cumsum_cone_size[j] + block_size[j];
    }
    if (cumsum_cone_size[k] != n)
        mexErrMsgTxt("block_size does not match X");

    if (m_in == 1)
    {
        // Allocate output matrix
        plhs[0] = mxCreateSparse(k, n, n, mxREAL);
        ir = mxGetIr(plhs[0]);
        jc = mxGetJc(plhs[0]);
        pr = mxGetPr(plhs[0]);
#ifndef SERIAL
#pragma omp parallel for private(j, i) shared(ir, jc, pr, X, k) schedule(dynamic, 1)
#endif
        for (j = 0; j < k; j++)
        {
            for (i = cumsum_cone_size[j]; i < cumsum_cone_size[j + 1]; i++)
            {
                ir[i] = j;
                jc[i] = i;
                pr[i] = X[i];
            }
        }
        jc[n] = n;
    }
    else
    {
        // Allocate output matrix
        plhs[0] = mxCreateSparse(n, k, n, mxREAL);
        ir = mxGetIr(plhs[0]);
        jc = mxGetJc(plhs[0]);
        pr = mxGetPr(plhs[0]);

        memcpy(jc, cumsum_cone_size, (k + 1) * sizeof(mwIndex));
#ifndef SERIAL
#pragma omp parallel for private(j, i) shared(ir, jc, pr, X, k) schedule(dynamic, 1)
#endif
        for (j = 0; j < k; j++)
        {
            for (i = cumsum_cone_size[j]; i < cumsum_cone_size[j + 1]; i++)
            {
                ir[i] = i;
                pr[i] = X[i];
            }
        }
    }

    mxFree(block_size);
    mxFree(cumsum_cone_size);
}
