%%  
%  Aalgo.AAt.At_col_denhor: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-08 10:17:15
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [d, algo] = invAAt(r, model, algo, params)

    % r can be matrix, each column is a r
    %% solve algo.AAt * d = r

    if params.precond_A
        d = r;
    elseif algo.AAt.isidentity
        d = r;
    else
        % _rd stands for reduced
        if model.n_box
            At_box = model.At_box;
            r_box = r(end - model.n_box + 1: end, :);
            r_rd = r(1: end - model.n_box, :) - 0.5 * At_box' * r_box;
        else
            r_rd = r;
        end
        
        if strcmp(algo.AAt.dense_column_strategy, 'augmented')
            [d] = solve_ldl([r_rd; zeros(algo.model.n_dense_col, size(r_rd, 2))], algo.AAt);
            d = d(1: end - algo.model.n_dense_col, :);
        elseif strcmp(algo.AAt.dense_column_strategy, 'SMW')
            % [d] = SMW(r_rd, algo.AAt.At_col_den, speye(size(algo.AAt.At_col_den, 1)), algo.AAt.inv_AAt);
            temp = algo.AAt.inv_AAt(r_rd);
            d = temp - algo.AAt.inv_AAt(algo.AAt.At_col_den' * (algo.AAt.smw_solve(algo.AAt.At_col_den * temp)));
        else %none
            [d] = algo.AAt.inv_AAt(r_rd);
        end



        if model.n_box
            d_box = 0.5 * (r_box - At_box * d);
            d = [d; d_box];
        end

        % fprintf("AAt residue = %e\n", norm(r -  AXfun(model.K, model.At_int, Atyfun(model.K, model.At_int, d))) / (1 + norm(r)));

    end

end