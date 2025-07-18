function [y] = SMW(x, Ut, invC, invD)
    %% use the Sherman-Morrison-Woodbury formula to compute the inverse of a matrix
    % (D + Ut' * C * Ut)^{-1} * x  = ( D^{-1} - D^{-1} * Ut' * (C + Ut * D^{-1} * Ut')^{-1} * Ut * D^{-1} ) * x
    % input:
    %   x: vector or matrix to be multiplied
    %   Ut: matrix 
    %   invC: inverse of C, matrix
    %   invD: inverse of D, functional handle
    if numel(Ut) == 0 || numel(invC) == 0
        y = invD(x);
        return;
    end
    UtinvDU = Ut * invD(Ut');
    UtinvDU = UtinvDU + invC;
    if nnz(UtinvDU) < 0.5 * numel(UtinvDU)
        UtinvDU = sparse(UtinvDU);
    else
        UtinvDU = full(UtinvDU);
    end
    invDx = invD(x);
    y = invDx - invD(Ut' * (UtinvDU \ (Ut * invDx)));
end
