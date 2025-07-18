function out = sparsity(obj)
    % sparsity of the whole cell 
    s = sum(cellfun(@numel, obj));
    if s == 0
        out = 0;
    else
        out = MatCell.nnz(obj) / s;
    end
end