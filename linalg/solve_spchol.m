function [d] = solve_spchol(rhs, info)
    assert(isfield(info, 'R'));
    assert(isfield(info, 'ordering'));
    d = info.R \ ( info.R' \ rhs(info.ordering) );
    d(info.ordering) = d;
end