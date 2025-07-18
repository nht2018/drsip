clear all; close all;
fprintf('-------------------------------------------------------------\n');
fprintf('Using Pardiso on a sparse, real and non-symmetric matrix.       \n');
fprintf('             Panua Technologies, www.panua.ch                \n');
fprintf('-------------------------------------------------------------\n');
% setenv("KMP_DUPLICATE_LIB_OK","TRUE")

% Computes the m solutions X to the collection of linear systems
%
%    A * X = B
%
% using the Pardiso solver, where A is a symmetric n x n matrix, B is an
% n x m matrix, and X is another n x m matrix.

% Specify parameters for the problem specification
n       = 100;
m       = n-1;
density = 0.7;

% Generate the matrix A 
A               = sprand(n,n,density) + rand(n,1).*speye(n);
[nrow_A,ncol_A] = size(A);
fprintf('Real non-symmetric matrix of size %d x %d, with %d nonzeroes \n',nrow_A,ncol_A,nnz(A));

% Generate a random collection of right-hand sides.
B = rand(n,m);

% Initialize the Pardiso internal data structures. Pardiso will
% handle the real non-symmetric matrices using the sparse direct solver.
% Arg1: 11 | Real and nonsymmetric matrix
% Arg2:  0 | Sparse direct solver
info = pardisoinit(11,0);

% Analyze the matrix and compute a symbolic factorization.
verbose = true;
info = pardisoreorder(A,info,verbose);
fprintf('The factors have %d nonzero entries.\n',info.iparm(18));

% Compute the numeric factorization.
info = pardisofactor(A,info,verbose);

% Compute the solutions X using the symbolic factorization.
[X, info] = pardisosolve(A,B,info,verbose);
fprintf('Pardiso performed %d iterative refinement steps.\n',info.iparm(7));

% Compute the residuals.
R = max(abs(A*X - B));
fprintf('The maximum residual for the solution X is %0.3g.\n',max(R(:)));

% Free the Pardiso data structures.
pardisofree(info);
clear info
