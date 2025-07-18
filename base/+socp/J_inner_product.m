function out = J_inner_product(X, Y, cone_size)
    % X0 * Y0 - bar X'* bar Y
    % to be implemented in mex in the future
    n = size(X, 1);
    ind_head = cumsum(cone_size) - cone_size + 1;
    out = 2 * X(ind_head) .* Y(ind_head) - blk_sum(X .* Y, cone_size, 1, false);
end