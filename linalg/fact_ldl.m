function [info, flag] = fact_ldl(mat, info, mat_regu)
    if nargin < 3 || isempty(mat_regu)
        mat_regu = speye(size(mat));
    end
    if nargin < 2 || isempty(info)
        info = struct();
    end
    if isfield(info, 'reg_init')
        beta = info.reg_init;
    else
        beta = 0;
    end
    % if condest(mat11) > 1e12 beta = 1e-12; end
    if beta
        mat = mat + beta * mat_regu;
    end
    [info.L, info.D] = ldl(mat);
    flag = 0;
end