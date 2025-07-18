%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-12 14:59:07
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-12 14:59:51
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [lhs] = linsysconcat(lhs, system_opt)
    lhs.mat22 = spdiag(lhs.mat22_diag);
    lhs.mat33 = spdiag(lhs.mat33_diag);

    lhs.dim1 = size(lhs.mat11, 1);
    lhs.dim2 = size(lhs.mat22, 1);
    lhs.dim3 = size(lhs.mat33, 1);

    if strcmp(system_opt, 'sparse') || strcmp(system_opt, 'sparse_psd')
        lhs.mat21 = [lhs.mat21_den; lhs.mat21_u];
    end

    lhs.mat = [lhs.mat11, lhs.mat21', lhs.mat31';
                lhs.mat21, lhs.mat22, lhs.mat32';
                lhs.mat31, lhs.mat32, lhs.mat33]; % acutually we need only lower part for factorization. this can be improved in the future

    lhs.lmut = @(x) lhs.mat * x;

    lhs.regu_diag = [ones(lhs.dim1, 1); sign(lhs.mat22_diag); sign(lhs.mat33_diag)];
    lhs.mat_regu = spdiag(lhs.regu_diag);
end