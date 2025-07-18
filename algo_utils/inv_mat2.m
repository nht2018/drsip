function [inv_func, info] = inv_mat2(mat, info)
% compute the inverse of a very sparse matrix.
% and return a function handle to solve the linear system. we split the solving step into two steps as forward and backward substitution
% After calling this function, the user can then call fwsolve2 and bwsolve2 to solve the linear system.
% we do not consider dense column handling here
% mat: a symmetric positive definite matrix
% inv_func: a functional handle, inverse of mat 
% info: a struct containing information about the factorization of mat


if nargin < 2
    info = struct();
end
if isdiag(mat)
    info.solver = 'isdiag';
    info.D = full(diag(mat));
    info.sqrtD = sqrt(full(diag(mat)));
    inv_func = @(x)  1 ./ info.D .* x;
else
    if nnz(mat) < 0.5 * numel(mat)
        info.solver = 'ldlchol';
        [info, flag] = fact_ldlchol(mat, info);
        [info.L, info.D] = ldlsplit(info.LD);
        info.sqrtD = sqrt(info.D);
        inv_func = @(x) solve_ldlchol(x, info);
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

    if flag ~= 0
        error("factorization failed in inv_mat");
    end
end


end