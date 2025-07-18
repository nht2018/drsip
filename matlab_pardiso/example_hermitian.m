clear all; close all;
fprintf('-------------------------------------------------------------\n');
fprintf('Using Pardiso on a complex Hermitian positive definite matrix.\n');
fprintf('             Panua Technologies, www.panua.ch                \n');
fprintf('-------------------------------------------------------------\n');

% Specify parameters for the problem specification
n       = 1e3;
lambda  = 3;

% Create the Hermitian positive definite matrix A and the vector b in the
% linear system Ax = b.
e = ones(n,1);
A = spdiags([ i*e lambda*e -i*e ],-1:1,n,n);
A = (A + A')/2;
b = randn(n,1);

[nrow_A,ncol_A] = size(A);
fprintf('Complex Hermitian matrix of size %d x %d, with %d nonzeroes \n',nrow_A,ncol_A,nnz(A));

% Compute solution to linear system
% ---------------------------------
% Initialize the Pardiso internal data structures. Pardiso will
% handle Hermitian positive definite matrices using the sparse direct solver. 
% Arg1: 4 | complex and Hermitian positive definite
% Arg2: 0 | sparse direct solver
info = pardisoinit(4,0);

% Analyze the matrix and compute a symbolic factorization.
verbose = true;
info = pardisoreorder(tril(A),info,verbose);
fprintf('The factors have %d nonzero entries.\n',info.iparm(18));

% Compute the numeric factorization.
info = pardisofactor(tril(A),info,verbose);

% Compute the solution v using the symbolic factorization.
[x, info] = pardisosolve(tril(A),b,info,true);
fprintf('Final solution computed.\n');
fprintf('Pardiso performed %d iterative refinement steps.\n',info.iparm(7));

% Compute the residuals.
r = abs(A*x - b);
fprintf('The maximum residual for the solution is %0.3g.\n',max(r));

% Free the Pardiso data structures.
pardisofree(info);

clear info;
