function [d] = solve_lu(rhs, info)
    assert(isfield(info, 'L'));
    assert(isfield(info, 'U'));
    d = info.U \ (info.L \ rhs);
end