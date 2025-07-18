function [inv_func, info] = inv_mat(Bt, coeff, D, info, dense_column_strategy)
    % compute the inverse of 
    % B * spdiag(coeff) * B' + D
    % here B donnot have many dense columns
    %      D is very sparse
    % input:
    %       Bt: the transpose of B
    %       coeff: a vector, the diagonal of the diagonal matrix
    %       dense_column_strategy: char, 'none', 'SMW', 'augmented'
    % usage:
    %       [inv_func, info] = inv_mat(Bt, coeff, D, info, dense_column_strategy)
    %       [inv_func, info] = inv_mat(Bt, 1, D, info, dense_column_strategy)
    %       [inv_func, info] = inv_mat(Bt, coeff, 0, info, dense_column_strategy)
    %       [inv_func, info] = inv_mat(Bt, coeff, D, [], dense_column_strategy)
    %       [inv_func, info] = inv_mat(Bt, coeff, D, info)
    %       [inv_func, info] = inv_mat(Bt, coeff, D)
    % inv_func: a functional handle, inverse of mat * spdiag(coeff) * mat' 
    % info: a struct containing information about the factorization of mat
    if nargin < 5
        dense_column_strategy = 'SMW';
    end
    if nargin < 4 || isempty(info)
        info = struct();
    end
    if nnz(Bt) < 0.5 * numel(Bt)
        Bt = sparse(Bt);
    else
        Bt = full(Bt);
    end

    if ~ issparse(Bt)
        dense_column_strategy = 'none' ;
    end
    assert(ismember(dense_column_strategy, {'none', 'SMW', 'augmented'}));
    if strcmp(dense_column_strategy, 'none')
        col_den = [];
        col_sp = [];
    else
        [col_den, col_sp] = detect_dense_row(Bt);     
    end


    if isempty(col_den) || isempty(col_sp) || strcmp(dense_column_strategy, 'none')
        if isscalar(coeff) && coeff == 1
            mat = Bt' * Bt;
        else
            mat = Bt' * spdiag(coeff) * Bt;
        end
        if ~ ( isscalar(D) && D == 0 )
            mat = mat + D;
        end
        if isdiag(mat)
            info.solver = 'isdiag';
            info.D = full(diag(mat));
            info.sqrtD = sqrt(full(diag(mat)));
            inv_func = @(x)  1 ./ info.D .* x;
        else
            if nnz(mat) < 0.5 * numel(mat)
                % info.solver = 'ldlchol';
                % [info, flag] = fact_ldlchol(mat, info);
                % inv_func = @(x) solve_ldlchol(x, info);
                info.solver = 'spchol';
                [info, flag] = fact_spchol(mat, info);
                inv_func = @(x) solve_spchol(x, info);
            else
                info.solver = 'chol';
                [info, flag] = fact_chol(mat, info);
                inv_func = @(x) solve_chol(x, info);
            end
            if flag ~=  0 % not positive definite
                if nnz(mat) < 0.5 * numel(mat)
                    info.solver = 'spldl';
                    [info, flag] = fact_spldl(mat, info);
                    inv_func = @(x) solve_spldl(x, info);
                else
                    info.solver = 'ldl';
                    [info, flag] = fact_ldl(mat, info);
                    inv_func = @(x) solve_ldl(x, info);
                end
            end
        end
    else
        col_den = intersect(col_den, find(coeff ~= 0));
        den = Bt(col_den, :);
        sp = Bt(col_sp, :);
        if isscalar(coeff) && coeff == 1
            square_sp = sp' * sp;  
        else
            square_sp = sp' * spdiag(coeff(col_sp)) * sp;
        end
        if ~ ( isscalar(D) && D == 0 )
            square_sp = square_sp + D;
        end
        if strcmp(dense_column_strategy, 'SMW')
            %% by SMW formula, 
            % mat^{-1} = (den' * spdiag(coeff_den) * den + square_sp)^{-1} 
            %                = square_sp^{-1} - square_sp^{-1} * den' * (spdiag(1 ./ coeff_den) + den * square_sp^{-1} * den')^{-1} * den * square_sp^{-1}
            if nnz(den) < 0.25 * numel(den)
                den = sparse(den);
            else
                den = full(den);
            end
            if isdiag(square_sp)
                info.solver = 'isdiag';
                diag_square_sp = full(diag(square_sp));
                assert(all(diag_square_sp >= 0));
                inv_square_sp = @(x) 1 ./ diag_square_sp .* x;
            else
                info.solver = 'ldlchol';
                [info, flag] = fact_ldlchol(square_sp, info);
                inv_square_sp = @(x) solve_ldlchol(x, info);
                if flag ~= 0
                    info.solver = 'spldl';
                    [info, flag] = fact_spldl(square_sp, info);
                    inv_square_sp = @(x) solve_spldl(x, info);
                end
            end
            % den_square_sp_inv_den is a small dense matrix since there is not many dense columns in At_sp
            inv_func = @(x) SMW(x, den, speye(size(den, 1)), inv_square_sp);
        else
            %% by augmented 
            % mat^{-1} rhs =  [squere_sp, den';            * [rhs;
            %                   den, diag(coeff_den)]^{-1}]    0]   
            den = sparse(den);
            % scaling = ones(numel(col_den), 1);
            scaling = sqrt(abs(coeff(col_den)));
            aug_mat = [square_sp, den' * spdiag(scaling); 
                      spdiag(scaling) * den, spdiag(- scaling.^ 2 ./ coeff(col_den))];
            info.solver = 'spldl';
            [info, flag] = fact_spldl(aug_mat, info);
            if flag ~= 0
                error('spldl factorization failed');
            end
            inv_func = @(x) extractFirstN(solve_spldl([x; zeros(numel(col_den), size(x, 2))], info), size(x, 1));
        end
    end

end

function x = extractFirstN(array, n)
    x = array(1:n, :);
end
