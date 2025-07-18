function [d] = solve_ldlchol(rhs, info)
    assert(isfield(info, 'ordering'));
    assert(isfield(info, 'LD'));
    d = ldlsolve(info.LD, rhs(info.ordering, :));
    d(info.ordering, :) = d;
end