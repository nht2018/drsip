function [d] = solve_splu(rhs, info)
    assert(isfield(info, 'L'));
    assert(isfield(info, 'U'));
    assert(isfield(info, 'P'));
    assert(isfield(info, 'Q'));
    d = zeros(size(rhs));
    d(info.Q, :) = info.U \ (info.L \ rhs(info.P, :));
end