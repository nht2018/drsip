%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-02-26 20:41:40
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-29 10:23:22
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [info, flag] = fact_pardiso(tril_mat, info, mat_regu)
    if nargin < 3 || isempty(mat_regu)
        mat_regu = speye(size(mat));
    end
    if nargin < 2 || isempty(info)
        info = struct();
    end
    if isfield(info, 'reg_init')
        beta = info.reg_init;
    else
        beta = 1e-16;
    end
    mat_temp = tril_mat + beta * mat_regu;
    if ~info.initialized
        info.pardiso_info = pardisoinit(-2, 0);
        info.initialized = true;
    end
    assert(isfield(info, 'pardiso_info'));
    info.pardiso_info = pardisoreorder(mat_temp, info.pardiso_info, false);
    % fprintf('[Pardiso] The factors have %d nonzero entries.\n',info.pardiso_info.iparm(18));
    info.pardiso_info = pardisofactor(mat_temp, info.pardiso_info, false);
    flag = 0;
    info.reg = beta;
end