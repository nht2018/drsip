function [d] = solve_ldl(rhs, info)
    assert(isfield(info, 'L'));
    assert(isfield(info, 'D'));
    assert(isfield(info, 'P'));
    assert(isfield(info, 'S'));
    d = info.S * rhs;
    d = d(info.P, :);
    d = info.L \ d; d = info.D \ d; d = info.L' \ d;
    d(info.P, :) = d;
    d = info.S * d;
end