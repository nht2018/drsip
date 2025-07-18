%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-11-26 20:37:29
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-26 20:38:00
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = diag_concat(A, B)
    out = [A, sparse(size(A, 1), size(B, 2));
             sparse(size(B, 1), size(A, 2)), B];
end