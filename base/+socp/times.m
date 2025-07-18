function out = times(X, Y, cone_size)
    % X \circ Y
    n = size(X, 1);
    ind_head = cumsum(cone_size) - cone_size + 1;
    out = repelem(X(ind_head), cone_size, 1) .* Y  + repelem(Y(ind_head), cone_size, 1) .* X;
    out(ind_head) = blk_sum(X .* Y, cone_size, 1, false);
end