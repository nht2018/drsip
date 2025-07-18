function [inv_mat22, schur_comp, info] = schur_complement(mat2t, coeff, mat11, mat12, info, dense_column_strategy)
    % compute the inverse of mat22 and schur complement of 
    % [mat11, mat12;
    % mat12', mat22]
    % where mat22 = mat2 * spdiag(coeff) * mat2' and mat2 is very sparse
    % note that size(mat2, 2) does not need to be the same as size(mat12, 2)
    %  coeff is a scalar or a vector of size size(mat2, 1), coeff need not to be positive
    % the input mat2t is the transpose of mat2
    % usage:
    %       [inv_mat22, schur_comp, info] = schur_complement(mat2t, coeff, mat11, mat12, dense_column_strategy)
    % inv_mat22: a functional handle, inverse of mat22
    % schur_comp = mat11 - mat12 * inv(mat22) * mat12'
    % info: a struct containing information about the factorization of mat22
    if nargin < 6
        dense_column_strategy = 'SMW';
    end

    [inv_mat22, info] = inv_mat(mat2t, coeff, 0, info, dense_column_strategy);
    
    schur_comp = mat11 - mat12 * inv_mat22(mat12');

end

