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
algo.AAt.dense_column_strategy = params.dense_column_strategy;
assert(ismember(algo.AAt.dense_column_strategy, {'default', 'none', 'SMW', 'augmented'}));
if ~ issparse(At_weighted) % the matrix is dense
    col_den = [];
    algo.AAt.dense_column_strategy = 'none';
    algo.AAt.isidentity = false;
else
    col_den = detect_dense_row(At_weighted);
    col_den = intersect(col_den, find(sum(At_weighted ~= 0, 2) > 1000));
    if isempty(col_den) && norm(At_weighted' * At_weighted - speye(size(At_weighted, 2)), 'fro') < 1e-12
        algo.AAt.isidentity = true;
    else
        algo.AAt.isidentity = false;
    end
    if params.precond_A % if we need to do preconditioning, we do not split the dense column
        col_den = [];
        algo.AAt.dense_column_strategy = 'none';
    end

end

if isempty(col_den) || numel(col_den) == size(At_weighted, 1) % no dense column or the whole matrix is dense
    algo.model.n_dense_col = 0;
    At0 = [];
    At1 = At_weighted;
else % otherwise we need to split dense column in computation
    fprintf(algo.fid, 'split %d dense columns, %d sparse columns\n', numel(col_den), size(At_weighted, 1) - numel(col_den));
    algo.model.n_dense_col = numel(col_den);
    At0 = At_weighted(col_den, :); % dense part
    At1 = At_weighted(setdiff(1:size(At_weighted, 1), col_den), :); % sparse part
end

%  detect very sparse rows in the matrix apart from the dense columns
if params.detect_very_sparse_rows
    row_sp = find(full(sum(At1 ~= 0, 1)) <= 2);
    if  numel(row_sp) > 10 && numel(row_sp) < size(At1, 2) - 10 % if the number of very sparse rows is small or almost all rows are very sparse number of very dense rows is large, we do not need to split the matrix
        algo.model.detect_very_sparse_rows = true;
    else
        algo.model.detect_very_sparse_rows = false;
    end
else
    algo.model.detect_very_sparse_rows = false;
end

if params.precond_A % Acutally we can do preconditioning with spliting dense and sparse rows. To be improved in the future
    algo.model.detect_very_sparse_rows = false;
end

%% decide the strategy to handle dense columns
if strcmp(algo.AAt.dense_column_strategy, 'default')
    if issparse(At_weighted)
        if numel(col_den) == 0 ||  numel(col_den) == size(At_weighted, 1) % no dense column 
            algo.AAt.dense_column_strategy = 'none' ;
        elseif numel(col_den) < 10 && ~ algo.model.detect_very_sparse_rows
            algo.AAt.dense_column_strategy = 'augmented' ;
        else
            algo.AAt.dense_column_strategy = 'SMW' ;
        end
    else
        algo.AAt.dense_column_strategy = 'none' ;
    end
end



assert(ismember(algo.AAt.dense_column_strategy, {'none', 'SMW', 'augmented'}) );
if strcmp(algo.AAt.dense_column_strategy, 'augmented')
    algo.model.detect_very_sparse_rows = false; % when we use strategy to handle dense columns, we do not handle very sparse rows at the same time; Actually we can do both, but it is not implemented yet.
    At0 = sparse(At0);
    AAt = [At1' * At1, At0';
        At0, -speye(size(At0, 1))];
    algo.AAt.solver = 'spldl';
    [algo.AAt, flag] = fact_spldl(AAt, algo.AAt);
    algo.AAt.time_fact = toc(t0) ;
    algo.AAt.dim = size(AAt, 1);
    algo.AAt.time_fact = toc(t0) ;
    
    return;
elseif strcmp(algo.AAt.dense_column_strategy, 'SMW')
    algo.model.detect_very_sparse_rows = false;
    % SMW
    % first we consider only the inverse of sparse part, i.e., (At1 * At1')^{-1}
    % and later use SMW formula to compute (At0' * At0 + At1' * At1)^{-1}
    if nnz(At0) < 0.25 * numel(At0)
        At0 = sparse(At0);
    else
        At0 = full(At0);
    end
    At1 = sparse(At1);
    
    
    algo.AAt.dim = size(At1, 2);
    empty_rows = find( sum(At1 ~= 0, 1) == 0 );
    nonempty_rows = setdiff(1: size(At1, 2), empty_rows);
    temp = sparse(numel(empty_rows), size(At0, 2));
    temp(:, empty_rows) = 1e-2 * speye(numel(empty_rows));
    temp_coeff0 = [ones(size(At0, 1), 1); - ones(numel(empty_rows), 1)];
    temp_coeff1 = [ones(size(At1, 1), 1); - ones(numel(empty_rows), 1)];
    At1 = [At1; temp];
    At0 = [At0; temp];
    if nnz(At0) < 0.25 * numel(At0)
        At0 = sparse(At0);
    else
        At0 = full(At0);
    end
    
    [inv_func, algo.AAt] = inv_mat(At1, temp_coeff1, 0, algo.AAt, 'none');
    algo.AAt.inv_AAt = inv_func;
    smw_temp = At0 * inv_func(At0');
    smw_temp = smw_temp + speye(size(smw_temp, 1));
    % [inv_func, algo.AAt] = inv_mat(At1, 1, 0, algo.AAt, 'none');
    % algo.AAt.inv_AAt = inv_func;
    % smw_temp = At0 * inv_func(At0');
    % smw_temp = smw_temp + speye(size(smw_temp, 1));
    % smw_temp = smw_temp + spdiag(temp_coeff0);
    
    if nnz(smw_temp) < 0.5 * numel(smw_temp)
        smw_temp = sparse(smw_temp);
        [smw_info, flag] = fact_ldlchol(smw_temp);
        algo.AAt.smw_solve = @(x) solve_ldlchol(x, smw_info);
    else
        smw_temp = full(smw_temp);
        [smw_info, flag] = fact_chol(smw_temp);
        algo.AAt.smw_solve = @(x) solve_chol(full(x), smw_info);
    end
    if flag ~= 0 % not positive definite
        if nnz(smw_temp) < 0.5 * numel(smw_temp)
            [smw_info, flag] = fact_spldl(smw_temp);
            algo.AAt.smw_solve = @(x) solve_spldl(x, smw_info);
        else
            [smw_info, flag] = fact_ldl(smw_temp);
            algo.AAt.smw_solve = @(x) solve_ldl(x, smw_info);
        end
    end
    if flag ~= 0
        error('SMW factorization failed');
    end
    algo.AAt.At_col_den = At0;
    algo.AAt.At_col_sp = At1;
    algo.AAt.smw_temp = smw_temp;
    algo.AAt.time_fact = toc(t0) ;
    return;
end

%% dense_column_strategy = 'none'
At1 = At_weighted; 
if algo.model.detect_very_sparse_rows
    row_den = setdiff(1:size(At1, 2), row_sp);
    At_sp = At1(:, row_sp);
    At_den = At1(:, row_den);
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
    
    [algo.AAt.inv_AAt_sp_sp, ~] = inv_mat(At_sp, 1, 0, struct(), 'SMW');
    if numel(col_den) > 0
        [inv_AAt_den_den, algo.AAt] = inv_mat(At_den_nonempty, 1, - AAt_den_sp * algo.AAt.inv_AAt_sp_sp(AAt_den_sp'), algo.AAt, 'SMW');
    else
        [inv_AAt_den_den, algo.AAt] = inv_mat(At_den_nonempty, 1, - AAt_den_sp * algo.AAt.inv_AAt_sp_sp(AAt_den_sp'), algo.AAt, 'none');
    end
    algo.AAt.inv_AAt_den_den = inv_AAt_den_den;
    algo.AAt.dim = size(At_den, 2);
    algo.AAt.inv_AAt = @(x) invAAt_temp(x, algo);
else
    algo.model.row_sp = [];
    algo.model.row_den = [1: size(At1, 2)];
    algo.model.detect_very_sparse_rows = false;
    algo.AAt.dim = size(At1, 2);
    [inv_func, algo.AAt] = inv_mat(At1, 1, 0, algo.AAt, 'none');
    algo.AAt.inv_AAt = inv_func;
end

algo.AAt.time_fact = toc(t0) ;
end


function [d] = invAAt_temp(r, algo)
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
