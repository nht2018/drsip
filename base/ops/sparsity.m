function out = sparsity(obj)
    % sparsity of the whole cell 
    if iscell(obj)
        s = sum(cellfun(@numel, obj));
        if s == 0
            out = 0;
        else
            out = nnz(obj) / s;
        end
    elseif isnumeric(obj)
        s = numel(obj);
        if s == 0
            out = 0;
        else
            out = nnz(obj) / s;
        end
    else
        error('input mUst be cell or matrix');
    end
end
