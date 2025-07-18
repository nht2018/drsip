%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-07 20:12:21
%  Description:
%
%  Copyright (c) 2023, Hantao Nie, Peking University.
%
function [model, algo] = factAAt(model, algo, params)
%% factorizes A * A' for the linear system A * A' * x = b
% the difference of factAAt2 lies in the strategy to handle very sparse rows

algo.AAt = struct;
t0 = tic;

%% ontain A in matrix form
At_weighted = model.At;
for p = 1: length(model.K)
    if strcmp(model.K{p}.type, 'b2l')
        At_weighted{p} = At_weighted{p} * sqrt(0.5);
    end
end

At_weighted = MatCell.vert_concat(At_weighted);
% split dense column and sparse column
if nnz(At_weighted) < 0.5 * numel(At_weighted)
    At_weighted = sparse(At_weighted);
else
    At_weighted = full(At_weighted);
end

%% detect dense column
% here we set the threshold to be 1000 because the multiplication and Cholesky factorization for dense matrix with size over 1000 work very slow on my computer, sometimes even cause memory error , hence we need to split the dense part and adress them carefully
algo.AAt.dense_column_strategy = 'none' ;
algo.AAt.isidentity = false;



%  detect very sparse rows in the matrix apart from the dense columns
if params.detect_very_sparse_rows
    row_sp = find(full(sum(At_weighted ~= 0, 1)) <= 2);
    if  numel(row_sp) > 10 && numel(row_sp) < size(At_weighted, 2) - 10 % if the number of very sparse rows is small or almost all rows are very sparse number of very dense rows is large, we do not need to split the matrix
        algo.model.detect_very_sparse_rows = true;
    else
        algo.model.detect_very_sparse_rows = false;
    end
else
    algo.model.detect_very_sparse_rows = false;
end





if algo.model.detect_very_sparse_rows
    row_den = setdiff(1:size(At_weighted, 2), row_sp);
    At_sp = At_weighted(:, row_sp);
    At_den = At_weighted(:, row_den);

    % AAt([row_den; row_sp], [row_den; row_sp]) = A_weighted * A_weighted'
    %     = [A_den * A_den'; A_den * A_sp';
    %        A_sp * A_den'; A_sp * A_sp']
    % note that A_sp * A_sp' is a diagonal matrix
    % solve the equation AAt * x = b is then equivalent to solve
    % (A_den * A_den' - A_den * A_sp' * (A_sp * A_sp')^{-1} * A_sp * A_den') * x_den = b_den - A_den * A_sp' * (A_sp * A_sp')^{-1} * b_sp
    % x_sp = (A_sp * A_sp')^{-1} ( b_sp - A_sp * A_den' x_den )
    % the matrix need to be factorized is A_den * A_den' - A_den * A_sp' * (A_sp * A_sp')^{-1} * A_sp * A_den'
    % and we need to store the A_den * A_sp' and A_sp * A_sp' for the calculation of x_sp
    nonempty_cols = find(sum(At_den ~= 0, 2) ~= 0);
    At_den_nonempty = At_den(nonempty_cols, :);
    if nnz(At_den_nonempty) < 0.25 * numel(At_den_nonempty)
        At_den_nonempty = sparse(At_den_nonempty);
    else
        At_den_nonempty = full(At_den_nonempty);
    end
    AAt_den_sp = At_den_nonempty' * At_sp(nonempty_cols, :);
    algo.AAt.AAt_den_sp = AAt_den_sp;
    algo.AAt.At_sp = At_sp;
    algo.AAt.AAt_den_den = At_den_nonempty' * At_den_nonempty ;
    algo.model.row_sp = row_sp;
    algo.model.row_den = row_den;
    algo.model.At_row_den = MatCell(length(model.At));
    algo.model.At_row_sp = MatCell(length(model.At));
    for p = 1: length(model.At)
        algo.model.At_row_den{p} = model.At{p}(:, algo.model.row_den);
        if nnz(algo.model.At_row_den{p}) < 0.25 * numel(algo.model.At_row_den{p})
            algo.model.At_row_den{p} = sparse(algo.model.At_row_den{p});
        else
            algo.model.At_row_den{p} = full(algo.model.At_row_den{p});
        end
        algo.model.At_row_sp{p} = model.At{p}(:, algo.model.row_sp);
    end
    
    AAt_sp_sp = At_sp' * At_sp;
    [algo.AAt.inv_AAt_sp_sp, algo.AAt.info_AAt_sp_sp] = inv_mat2(AAt_sp_sp, struct);

    [inv_AAt_den_den, algo.AAt] = inv_mat2(At_den_nonempty' * At_den_nonempty - AAt_den_sp * algo.AAt.inv_AAt_sp_sp(AAt_den_sp'), algo.AAt);

    algo.AAt.inv_AAt_den_den = inv_AAt_den_den;
    algo.AAt.dim = size(At_den, 2);
    algo.AAt.inv_AAt = @(x) invAAt_2blocks(x, algo);
    algo.AAt.fwsolve = @(x) fwsolve_2blocks(x, algo);
    algo.AAt.bwsolve = @(x) bwsolve_2blocks(x, algo);
else
    algo.model.row_sp = [];
    algo.model.row_den = [1: size(At_weighted, 2)];
    algo.model.detect_very_sparse_rows = false;
    algo.AAt.dim = size(At_weighted, 2);
    [inv_func, algo.AAt] = inv_mat(At_weighted, 1, 0, algo.AAt, 'none');
    if strcmp(algo.AAt.solver, 'ldlchol')
        [algo.AAt.L, algo.AAt.D] = ldlsplit(algo.AAt.LD);
        algo.AAt.D = diag(algo.AAt.D);
        algo.AAt.sqrtD = sqrt(algo.AAt.D);
    end
    algo.AAt.inv_AAt = inv_func;
    algo.AAt.fwsolve = @(x) fwsolve2(x, algo.AAt);
    algo.AAt.bwsolve = @(x) bwsolve2(x, algo.AAt);
end

algo.AAt.time_fact = toc(t0) ;
end


function [d] = invAAt_2blocks(r, algo)
    row_sp = algo.model.row_sp;
    row_den = algo.model.row_den;
    AAt_den_sp = algo.AAt.AAt_den_sp;
    r_sp = r(row_sp, :);
    r_rd = r(row_den, :) - AAt_den_sp * algo.AAt.inv_AAt_sp_sp(r_sp );
    [d_den] = algo.AAt.inv_AAt_den_den(r_rd) ;
    d_sp = algo.AAt.inv_AAt_sp_sp(r_sp - AAt_den_sp' * d_den);
    d = zeros(size(r_rd));
    d(row_den, :) = d_den;
    d(row_sp, :) = d_sp;
end


function [d] = fwsolve_2blocks(r, algo)

    row_sp = algo.model.row_sp;
    row_den = algo.model.row_den;
    AAt_den_sp = algo.AAt.AAt_den_sp;
    r_den = r(row_den, :);
    r_sp = r(row_sp, :);
    r_den = r_den - AAt_den_sp * algo.AAt.inv_AAt_sp_sp(r_sp);
    [d_den] = fwsolve2(r_den, algo.AAt);
    [d_sp] = fwsolve2(r_sp, algo.AAt.info_AAt_sp_sp);
    d = zeros(size(r));
    d(row_den, :) = d_den;
    d(row_sp, :) = d_sp;
end


function [d] = bwsolve_2blocks(r, algo)
    row_sp = algo.model.row_sp;
    row_den = algo.model.row_den;
    AAt_den_sp = algo.AAt.AAt_den_sp;
    r_den = r(row_den, :);
    r_sp = r(row_sp, :);
    [d_den] = bwsolve2(r_den, algo.AAt);
    [d_sp] = bwsolve2(r_sp, algo.AAt.info_AAt_sp_sp);
    d_sp = d_sp - algo.AAt.inv_AAt_sp_sp(AAt_den_sp' * d_den);
    d = zeros(size(r));
    d(row_den, :) = d_den;
    d(row_sp, :) = d_sp;
end
