%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-27 12:46:57
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [lhs] = schurmat(lhs, H, p, model, algo)
    %% system_opt
    % 1: ldl on sparse matrix
    % 3: chol on dense matrix 
    % 4: iterative method
    At = model.At;
    system_opt = algo.newton.system_opt;
    
    if system_opt == 1
        At_den = algo.model.At_col_den{p};
        At_sp = algo.model.At_col_sp{p};
        idxden = algo.model.col_den{p};
        idxsp = algo.model.col_sp{p};
        lhs.mat11 = lhs.mat11 + At_sp' * spdiag( H{p}.shift(idxsp)) * At_sp;
        lhs.mat31 = [lhs.mat31; At_den];
        lhs.mat33 = diag_concat(lhs.mat33, spdiag(- 1 ./ H{p}.shift(idxden))) ;
    elseif system_opt == 2
        At_den = algo.model.At_row_den{p};
        At_sp = algo.model.At_row_sp{p};
        nonempty_cols = find(sum(At_den ~= 0, 2) ~= 0);
        At_den_nonempty = At_den(nonempty_cols, :);
        if nnz(At_den_nonempty) < 0.5 * numel(At_den_nonempty)
            At_den_nonempty = sparse(At_den_nonempty);
        else
            At_den_nonempty = full(At_den_nonempty);
        end
        lhs.mat11 = lhs.mat11 + At_den_nonempty' * spdiag(H{p}.shift(nonempty_cols)) * At_den_nonempty;
        lhs.mat12 = lhs.mat12 + At_den_nonempty' * spdiag(H{p}.shift(nonempty_cols)) * At_sp(nonempty_cols, :);
        % lhs.mat22 = lhs.mat22 + At_sp' * spdiag(H{p}.shift) * At_sp;
        lhs.mat2t = [lhs.mat2t; At_sp];
        lhs.mat22_coeff = [lhs.mat22_coeff; H{p}.shift];
    elseif system_opt == 3
        lhs.mat11 = lhs.mat11 + At{p}' * spdiag(H{p}.shift) * At{p};
    elseif system_opt == 4
        lhs.mat11{p} = @(r_) At{p}' * (H{p}.shift .* (At{p} * r_) );
    end

end