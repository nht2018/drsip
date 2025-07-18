function [d] = solve_pardiso(lhs, rhs, info)
    assert(isfield(info, 'pardiso_info'));
    [d, info.pardiso_info] = pardisosolve(lhs.tril_mat, rhs, info.pardiso_info, false);
end