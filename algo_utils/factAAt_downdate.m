%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-31 17:54:02
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [model, algo] = factAAt_downdate(downdate_mat, model, algo, params)
    %% factorizes A * A' for the linear system A * A' * x = b

    if isfield(model, 'AAt')
        AAt = model.AAt;
    else
        AAt = AAtfun(model.At) ;
    end
    if algo.AAt.isidentity;
        error('to be implemented in the future');
    end
    if strcmp(algo.AAt.solver, 'ldlchol')
        [algo.AAt.LD] = ldlupdate(algo.AAt.LD, downdate_mat(algo.AAt.perm, :), '-');
        [algo.AAt.L, algo.AAt.D] = ldlsplit(algo.AAt.LD);
        algo.AAt.Lt = algo.AAt.L';
        algo.AAt.D = sqrt(algo.AAt.D);
    elseif strcmp(algo.AAt.solver, 'lchol')
        error('lchol do not support downdate');
    elseif strcmp(algo.AAt.solver, 'spchol')
        error('spchol do not support downdate');
    elseif strcmp(algo.AAt.solver, 'chol')
        [algo.AAt.R, flag] = cholupdate(algo.AAt.R, downdate_mat, '-');
        if flag ~= 0
            error('cholupdate failed');
        end
        algo.AAt.Rt = algo.AAt.R';
    else
        error('unknown solver for AAt');
    end
    
end