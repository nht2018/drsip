% pardisoinit.m Help file for the pardisoinit MEX-file.
%  Interface between the data initialization in the Panua-Pardiso 8.0 solver and Matlab.
% 
%  The calling syntax is:
% 
%  >> info = pardisoinit(Matrix_Type, Solver);
% 
% Matrix_Type:
%    1  real and structurally symmetric
%    2  real and symmetric positive definite
%   -2  real and symmetric indefinite
%    3  complex and structurally symmetric
%    4  complex and Hermitian positive definite
%   -4  complex and Hermitian indefinite 6 complex and symmetric
%   11  real and nonsymmetric
%   13  complex and nonsymmetric
% 
% Solver:
%    0  sparse direct solver
%    1  multi-recursive iterative solver
% 
% For more information, refer to the Panua-Pardiso manual, available at: https://panua.ch/pardiso.
% 
%  Copyright (C) 2022 until present, by Panua Technologies Sagl, Switzerland. All Rights Reserved.                       
%  This program can be downloaded from: https://panua.ch. 
