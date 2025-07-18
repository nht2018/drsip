%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-02-26 20:41:40
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-29 10:22:52
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [info, flag] = fact_chol(mat, info, mat_regu)
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
    n_iter = 0;
    while true
        mat_temp = mat + beta * mat_regu;
        [info.R, flag] = chol(mat_temp);
        n_iter = n_iter + 1;
        if flag == 0 || n_iter > 10
            break;
        elseif beta == 0
            beta = 1e-12;
            % beta = 1e-16 * norm(model.Atbar, 'fro');
        else
            beta = beta * 10;
        end
    end
    info.reg = beta;
end
