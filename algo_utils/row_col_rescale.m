%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-19 11:49:43
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [model, algo] = row_col_rescale(model, algo, option)
    %% rescale the data matrix A so that its norm approximate 1
    % Input
    %   model: generate socp model
    %   option = 1: 1 norm
    %               T. Pock and A. Chambolle. Diagonal preconditioning for first order primal-dual algorithms in
    %               convex optimization. In 2011 International Conference on Computer Vision, pages 1762–1769.
    %               IEEE, 2011.                        
    %            2:  2 norm
    %            3:  infty norm
    %               D. Ruiz. A scaling algorithm to equilibrate both rows and columns norms in matrices. Technical report, CM-P00040415, 2001
    % 
    % Output
    %  model: rescaled model
    %  algo.rescale.left_scale_A: a vector of size  [m, 1]
    %  algo.rescale.left_scale_F: a vector of size  [k, 1]
    %                       for rows that in the same secong order cone block, they have the same scaling factor
    %  algo.rescale.right_scale: a vector of size [n, 1]
    %
    % explanations on some intermediate variables
    % row_norm :  store the norm of each row in  [A; F], size = [m + k, 1]
    % col_norm :  store the norm of each column in [A; F], size = [n, 1]

    t0 = tic;
    if nargin < 2
        option = 2;
    end
    % fprintf("[before rescale] A norm: %e F norm: %e \n", norm(model.A, 'fro'), norm(model.F, 'fro'));
    algo.rescale = struct; 
    if ~isfield(model, 'At_mat')
        At = MatCell.vert_concat(model.At);
    else
        At = model.At_mat;   
    end
    [n, m] = size(At);
    left_scale = ones(m, 1);
    right_scale = ones(n, 1);


    n_iter = 5; % number of iterations. For large model, each iteration is time consuming. Hence we only use a few iteration .


    K = model.K;

    % idx_cone = []; % we dont want to scale columns of conic variables
    % ind = 1;
    % for p = 1: length(K)
    %     cone = K{p};
    %     if strcmp(cone.type, 'l') || strcmp(cone.type, 'r') || strcmp(cone.type, 's')
    %         idx_cone = [idx_cone, ind: ind + sum(cone.size) - 1];
    %     end
    %     ind = ind + sum(cone.size);
    % end

    for i = 1: n_iter
        if option == 1  % use 1 norm
            [col_norm, row_norm] = row_col_norm(At, 1, 1);
            col_norm = conic_average(col_norm, K);
        elseif option == 2 % use 2 norm
            [col_norm, row_norm] = row_col_norm(At, 2, 2);
            col_norm = conic_average(col_norm, K);
        elseif option == 3 % use infty norm
            [col_norm, row_norm] = row_col_norm(At, 3, 3);
            col_norm = conic_average(col_norm, K);
        end
        temp_left = 1 ./ sqrt(row_norm);
        temp_right = 1 ./ sqrt(col_norm);
        At = spdiag(temp_right) * At * spdiag(temp_left);

        left_scale = left_scale .* temp_left;
        right_scale = right_scale .* temp_right;
    end

    blk_size = cellfun(@(x) size(x, 1), model.At);
    algo.rescale.left_scale = left_scale;
    algo.rescale.right_scale = mysmat(K, MatCell.vert_split(right_scale, blk_size));
    model.At_mat = At;
    model.At = MatCell.vert_split(model.At_mat,  blk_size);
    model.b = model.b .* algo.rescale.left_scale;
    % x = 1./ right_scale .* x;
    model.c = model.c .* algo.rescale.right_scale;
    for p = 1: length(K)
        cone = K{p};
        if strcmp(cone.type, 'b2l') || strcmp(cone.type, 'b')
            model.K{p}.params = model.K{p}.params ./ algo.rescale.right_scale{p};
        end
    end

    % fprintf("[after rescale]  A norm: %e F norm: %e \n", norm(model.A, 'fro'), norm(model.F, 'fro'));
    algo.rescale.time = toc(t0);

end


function [col_norm_] = conic_average(col_norm_, K)
    ind = 1;
    for p = 1: length(K)
        cone = K{p};
        if strcmp(cone.type, 'q') || endsWith(cone.type, 'q') || strcmp(cone.type, 'r') || strcmp(cone.type, 's')
            col_norm_(ind: ind + sum(cone.size) - 1) = blk_average(col_norm_(ind: ind + sum(cone.size) - 1), cone.size, 1, true);
        end
        ind = ind + sum(cone.size);
    end
end