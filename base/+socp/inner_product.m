function out = inner_product(X, Y, cone_size)
    % X dot Y
    % to be implemented in mex in the future
    out = blk_sum(X .* Y, cone_size, 1, false);
end