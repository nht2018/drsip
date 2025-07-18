%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-02-26 20:41:40
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-05 15:52:01
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   

function [info, flag] = fact_spldl(mat, info, mat_regu)
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
    % if condest(mat11) > 1e12 beta = 1e-12; end
    if beta
        mat = mat + beta * mat_regu;
    end
    [info.L, info.D, info.P, info.S] = ldl(mat, 'vector');
    flag = 0;
end