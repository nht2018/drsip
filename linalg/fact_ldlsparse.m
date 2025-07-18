%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-02-26 20:41:40
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-29 10:23:08
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [info, flag] = fact_ldlsparse(mat, info, mat_regu)
    if nargin < 3 || isempty(mat_regu)
        mat_regu = speye(size(mat));
    end
    if nargin < 2 || isempty(info)
        info = struct();
    end
    if isfield(info, 'reg_init')
        beta = info.reg_init;
    else
        beta = 0;
    end
    mat = sparse(mat);
    mat_temp = mat + beta * mat_regu;
    if ~isfield(info, 'ordering') || isempty(info.ordering)
        info.ordering = analyze(mat_temp, 'sym');
    end
    [info.L, info.D, parent, flops] = ldlsparse(mat_temp, info.ordering);
    info.L = info.L + speye(size(info.L));
    info.reg = beta;
    flag = 0; % to be improved
end