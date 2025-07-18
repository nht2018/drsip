
function [info, flag] = fact_splu(mat, info, mat_regu)
    if nargin < 3 || isempty(mat_regu)
        mat_regu = speye(size(mat));
    end
    if nargin < 2 || isempty(info)
        info = struct();
    end
    if ~ issparse(mat)
        mat = sparse(mat);
    end
    if isfield(info, 'reg_init')
        beta = info.reg_init;
    else
        beta = 0;
    end
    [info.L, info.U, info.P, info.Q] = lu(mat + beta * mat_regu, 'vector');
    flag = 0;
end