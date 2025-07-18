function [d] = solve_lchol(rhs, info)
    assert(isfield(info, 'L'));
    d = info.L' \ ( info.L \ rhs );
end