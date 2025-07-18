function [d] = solve_ldlsparse(rhs, info)
    assert(isfield(info, 'L'));
    assert(isfield(info, 'D'));
    assert(isfield(info, 'ordering'));
    d = info.L' \ (info.D \ (info.L \ rhs(info.ordering, :)));
    d(info.ordering, :) = d;
end