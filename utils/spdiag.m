function out = spdiag(input)
    if length(input) == 0
        out = sparse(0,0);
    else
        assert(size(input, 1) == 1 || size(input, 2) == 1)
        out = spdiags(input(:),0,length(input),length(input));
    end
    
end