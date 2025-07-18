function C = unaryOperation(A, opname)
    if iscell(A)
        C = cellfun(@(a) builtin(opname, a), A, 'UniformOutput', false);
    else
        C = builtin(opname, A);
    end
end