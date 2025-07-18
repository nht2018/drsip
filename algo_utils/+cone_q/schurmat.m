%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-27 12:45:09
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
    cone = model.K{p};
    system_opt = algo.newton.system_opt;
    if system_opt == 1
        At_den = algo.model.At_col_den{p};
        At_sp = algo.model.At_col_sp{p};
    elseif system_opt == 2
        At_den = algo.model.At_row_den{p};
        At_sp = algo.model.At_row_sp{p};
    end

    if size(H{p}.shift, 1) == length(cone.size)
        H{p}.shift = repelem(H{p}.shift, cone.size, 1); % size = [n, 1]
    end

    if size(H{p}.lr, 2) < size(H{p}.coeff, 1)
        Hlr_cell = arrayfun(@(j)  blk_spdiag(H{p}.lr(:, j), cone.size), 1:size(H{p}.lr, 2), 'UniformOutput', false);
        Hlr = horzcat(Hlr_cell{:}) ;
    else
        Hlr = H{p}.lr ;
    end
    if system_opt ~= 2
        AHlr = At{p}' * Hlr ;
    else
        AdenHlr = At_den' * Hlr ;
        AspHlr = At_sp' * Hlr ;
    end

    if size(H{p}.coeff, 2) == 1 
        coeff_zeros = find(H{p}.coeff == 0);
    end
    
    if system_opt == 1
        % low rank part
        % denote by [m, n] = size(A), k = length(cone.size)
        
        [idxden, idxsp] = detect_dense_col(AHlr);
        if exist('coeff_zeros', 'var')
            % idxsp = setdiff(idxsp, coeff_zeros);
            idxsp = union(idxsp, coeff_zeros);
            idxden = setdiff(idxden, coeff_zeros);
        end

        if ~ isempty(idxden) % separate sparse and dense handling
            if ~ isempty(idxsp)
                AHlr_sp = AHlr(:, idxsp);
                if size(H{p}.coeff, 2) == 1 
                    lhs.mat11 = lhs.mat11 + AHlr_sp * spdiag(H{p}.coeff(idxsp)) * AHlr_sp';
                else
                    lhs.mat11 = lhs.mat11 + AHlr_sp * H{p}.coeff(idxsp, idxsp) * AHlr_sp';
                end
            end

            if size(H{p}.coeff, 2) == 1 
                lhs.mat31 = [lhs.mat31; AHlr(:, idxden)'];
                if isfield(lhs, 'mat33_diag')
                    lhs.mat33_diag = [lhs.mat33_diag; - 1 ./ H{p}.coeff(idxden)];
                else
                    lhs.mat33 = diag_concat(lhs.mat33, spdiag(- 1 ./ H{p}.coeff(idxden)));         
                end
            else
                lhs.mat31 = [lhs.mat31; AHlr(:, idxden)'];
                if ~ isfield(H{p}, 'invcoeff')
                    H{p}.invcoeff = inv(H{p}.coeff);
                end
                assert( ~ isfield(lhs, 'mat33_diag') )
                lhs.mat33 = diag_concat(lhs.mat33, - H{p}.invcoeff(idxden, idxden)); 
            end
        else % all columns are sparse
            if size(H{p}.coeff, 2) == 1 
                Hpcoeff = spdiag(H{p}.coeff);
            else
                Hpcoeff = H{p}.coeff;
            end
            lhs.mat11 = lhs.mat11 + AHlr * Hpcoeff * AHlr';
        end

        % diagonal part
        
        idxden = algo.model.col_den{p};
        idxsp = algo.model.col_sp{p};
    
        lhs.mat11 = lhs.mat11 + At_sp' * spdiag( H{p}.shift(idxsp)) * At_sp;
        if ~isempty(idxden)
            lhs.mat31 = [lhs.mat31; At_den];
            if isfield(lhs, 'mat33_diag')
                lhs.mat33_diag = [lhs.mat33_diag; - 1 ./ H{p}.shift(idxden)];
            else
                lhs.mat33 = diag_concat(lhs.mat33, spdiag(- 1 ./ H{p}.shift(idxden)));
            end
        end
    elseif system_opt == 2
        nonempty_cols = find(sum(At_den ~= 0, 2) ~= 0);
        At_den_nonempty = At_den(nonempty_cols, :);
        if nnz(At_den_nonempty) < 0.5 * numel(At_den_nonempty)
            At_den_nonempty = sparse(At_den_nonempty);
        else
            At_den_nonempty = full(At_den_nonempty);
        end
        if size(H{p}.coeff, 2) == 1 
            Hpcoeff = spdiag(H{p}.coeff);
        else
            error("to be implemented for matrix form H{p}.coeff")
        end
        lhs.mat11 = lhs.mat11 + AdenHlr * Hpcoeff * AdenHlr' + At_den_nonempty' * spdiag(H{p}.shift(nonempty_cols)) * At_den_nonempty;
        lhs.mat12 = lhs.mat12 + AdenHlr * ( Hpcoeff * AspHlr') + At_den_nonempty' * spdiag(H{p}.shift(nonempty_cols)) * At_sp(nonempty_cols, :);
        % lhs.mat22 = lhs.mat22 + AspHlr * Hpcoeff * AspHlr' + At_sp' * spdiag(H{p}.shift) * At_sp;
        lhs.mat2t = [lhs.mat2t; AspHlr';  At_sp];
        lhs.mat22_coeff = [lhs.mat22_coeff; H{p}.coeff; H{p}.shift];
    elseif  system_opt == 3
        % low rank part
        if size(H{p}.coeff, 2) == 1 
            Hpcoeff = spdiag(H{p}.coeff);
        else
            Hpcoeff = H{p}.coeff;
        end
        lhs.mat11 = lhs.mat11 + AHlr * Hpcoeff * AHlr';

        % diagonal part
        lhs.mat11 = lhs.mat11 + At{p}' * spdiag(H{p}.shift) * At{p};
    elseif system_opt == 4
        % low rank part
        if size(H{p}.coeff, 2) == 1 
            lhs.mat11{p} = @(r_) At{p}' * ( spdiag(H{p}.shift) * (At{p} * r_) )  + AHlr * (spdiag(H{p}.coeff) * (AHlr' * r_) );
        else
            lhs.mat11{p} = @(r_) At{p}' * ( spdiag(H{p}.shift) * (At{p} * r_) )  + AHlr * (H{p}.coeff * (AHlr' * r_) );
        end
    end

end