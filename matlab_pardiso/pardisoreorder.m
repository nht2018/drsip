% pardisoreorder.m Help file for the pardisoreorder MEX-file.
%  Interface between the reordering techniques in the Panua-Pardiso 8.0 solver and Matlab.
% 
%  The calling syntax is:
% 
%  >> info = pardisoreorder(tril(A), info, verbose);
% 
%  Users can supply their own (vector) reordering p by:
% 
%  >> info = pardisoreorder(tril(A), info, verbose, p);
% 
%  For more information, refer to the Panua-Pardiso manual, available at: https://panua.ch/pardiso.
% 
%  Copyright (C) 2022 until present, by Panua Technologies Sagl, Switzerland. All Rights Reserved.                       
%  This program can be downloaded from: https://panua.ch. 
