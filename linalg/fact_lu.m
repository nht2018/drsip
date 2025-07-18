%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-02-26 20:41:40
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-29 10:23:15
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [info, flag] = fact_lu(mat, info, mat_regu)
    if nargin < 3 || isempty(mat_regu)
        mat_regu = speye(size(mat));
    end
    if nargin < 2 || isempty(info)
        info = struct();
    end
    if issparse(mat)
        mat = full(mat);
    end
    if isfield(info, 'reg_init')
        beta = info.reg_init;
    else
        beta = 0;
    end
    [info.L, info.U] = lu(mat + beta * mat_regu);
    flag = 0;
end