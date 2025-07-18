%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-03-01 14:42:16
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-01 14:42:20
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [inv_func, schur_comp] = schur_complement(mat2t, coeff, mat11, mat12, dense_column_strategy)
    % compute the schur complement of 
    % [mat11, mat12;
    % mat12', mat22]
    % where mat22 = mat2 * spdiag(coeff) * mat2' and mat2 is very sparse
    % note that size(mat2, 2) does not need to be the same as size(mat12, 2)
    %  coeff is a scalar or a vector of size size(mat2, 1), coeff need not to be positive
    % the input mat2t is the transpose of mat2
    % [inv_func, schur_comp] = schur_complement(mat11, mat12, mat2)
    % inv_func: a functional handle, inverse of mat22
    % schur_comp = mat11 - mat12 * inv(mat22) * mat12'
    if nargin < 5
        dense_column_strategy = 'SMW';
    end
    assert(ismember(dense_column_strategy, {'none', 'SMW', 'augmented'}));
    [col_den, col_sp] = detect_dense_row(mat2t);

    if isscalar(coeff)
        coeff = coeff * ones(size(mat2t, 1), 1);
    end

    if isempty(col_den) || strcmp(dense_column_strategy, 'none')
        mat22 = mat2t' * spdiag(coeff) * mat2t;
        if isdiag(mat22)
            diag_mat22 = full(diag(mat22));
            assert(all(diag_mat22 >= 0));
            schur_comp = mat11 - mat12 * spdiag(1 ./ diag_mat22) * mat12';
            inv_func = @(x) 1 ./ diag_mat22 .* x;
        else
            L = chol(mat22, 'lower');
            temp = L \ mat12';
            schur_comp = mat11 - temp' * temp;
            inv_func = @(x) L' \ (L \ x);
        end
    else
        den = mat2t(col_den, :);
        sp = mat2t(col_sp, :);
        square_sp = sp' * spdiag(coeff(col_sp)) * sp;    
        if strcmp(dense_column_strategy, 'SMW')
            %% by SMW formula, 
            % mat22^{-1} = (den' * spdiag(coeff_den) * den + square_sp)^{-1} 
            %                = square_sp^{-1} - square_sp^{-1} * den' * (spdiag(1 ./ coeff_den) + den * square_sp^{-1} * den')^{-1} * den * square_sp^{-1}
            if nnz(den) < 0.5 * numel(den)
                den = sparse(den);
            else
                den = full(den);
            end
            if isdiag(square_sp)
                diag_square_sp = full(diag(square_sp));
                assert(all(diag_square_sp >= 0));
                inv_square_sp = @(x) 1 ./ diag_square_sp .* x;
                den_square_sp_inv_den = den * spdiag(1 ./ diag_square_sp) * den';
            else
                L = chol(square_sp, 'lower');
                inv_square_sp = @(x) L' \ (L \ x);
                temp = L \ den';
                den_square_sp_inv_den = temp' * temp;
            end
            % den_square_sp_inv_den is a small dense matrix since there is not many dense columns in At_sp
            den_square_sp_inv_den = full(den_square_sp_inv_den) + spdiag(1 ./ coeff(col_den));

            inv_func = @(x) inv_square_sp(x) - inv_square_sp(den' * ( den_square_sp_inv_den \ (den * inv_square_sp(x))));
        else
            %% by augmented 
            % mat22^{-1} rhs =  [squere_sp, den';            * [rhs;
            %                   den, diag(coeff_den)]^{-1}]    0]   
            den = sparse(den);
            % scaling = ones(numel(col_den), 1);
            scaling = sqrt(abs(coeff(col_den)));
            aug_mat = [square_sp, den' * spdiag(scaling); 
                      spdiag(scaling) * den, spdiag(- scaling.^ 2 ./ coeff(col_den))];
            info = struct;
            [info, flag] = fact_ldl(aug_mat, info);
            if flag ~= 0
                error('LDL factorization failed');
            end
            inv_func = @(x) extractFirstN(solve_ldl([x; zeros(numel(col_den), size(x, 2))], info), size(x, 1));
        end
        schur_comp = mat11 - mat12 * inv_func(mat12');
    end

end


function x = extractFirstN(array, n)
    x = array(1:n, :);
end