function C = elementwiseOperation(A, B, opname)
    if iscell(A) && iscell(B)
        C = cellfun(@(a, b) builtin(opname, a, b), A, B, 'UniformOutput', false);
    elseif iscell(A) 
        C = cellfun(@(a) builtin(opname, a, B), A, 'UniformOutput', false);
    elseif iscell(B)
        C = cellfun(@(b) builtin(opname, A, b), B, 'UniformOutput', false);
    else
        C = builtin(opname, A, B);
    end
end