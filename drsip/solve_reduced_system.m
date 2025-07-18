%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-02-25 17:26:31
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-07 21:09:31
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [d, res_rd, algo] = solve_reduced_system(lhs, rhs, model, algo, params)
    if algo.newton.system_opt == 1
        % d = lhs.mat  \ rhs;
        if strcmp(algo.newton.dense_column_strategy, 'augmented')
            [d, algo.newton] = linsyssolve(lhs, rhs, algo.newton);
            res_rd = norm(rhs - lhs.lmut(d)) / (1 + norm(rhs)) ;
            d = d(1: lhs.dim1 + lhs.dim2, :);
        elseif strcmp(algo.newton.dense_column_strategy, 'SMW')
            invAKAt = @(x) linsyssolve(lhs, x, algo.newton);
            d = SMW(rhs, [lhs.mat21, lhs.mat31], - lhs.mat33, invAKAt);
            res_rd = norm(rhs - lhs.lmut(d)) / (1 + norm(rhs)) ;
        elseif strcmp(algo.newton.dense_column_strategy, 'none') % same to system_opt = 3
            [d, algo.newton] = linsyssolve(lhs, rhs, algo.newton);
            res_rd = norm(rhs - lhs.lmut(d)) / (1 + norm(rhs)) ;
        else
            error('Unknown dense column strategy %s', algo.newton.dense_column_strategy);
        end
    elseif algo.newton.system_opt == 2
        row_sp = algo.model.row_sp;
        row_den = algo.model.row_den;
        rhs_sp = rhs(row_sp, :);
        rhs_rd = rhs(row_den, :) - lhs.mat12 * lhs.inv_mat22(rhs_sp);
        [d_den, algo.newton] = linsyssolve(lhs, rhs_rd, algo.newton);
        d_sp = lhs.inv_mat22(rhs_sp - lhs.mat12' * d_den);
        d = zeros(size(rhs));
        d(row_den, :) = d_den;
        d(row_sp, :) = d_sp;
        res_rd = norm(rhs - lhs.lmut(d)) / (1 + norm(rhs)) ;
    elseif ismember(algo.newton.system_opt, [3, 4])
        % d = lhs.mat  \ rhs;
        [d, algo.newton] = linsyssolve(lhs, rhs, algo.newton);
        res_rd = norm(rhs - lhs.lmut(d)) / (1 + norm(rhs)) ;
    else
        error('Unknown system option');
    end

    if any(isnan(d))
        warning('d contains NaN');
    end
%% check residual of reduced system



end