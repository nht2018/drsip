%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-19 11:51:13
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [At, b, c, lb, ub, rescale] = row_col_rescale_mat(At, b, c, lb, ub, option)
    %% rescale the data matrix A so that its norm approximate 1.
    % Input
    %   At: a matrix of size [n, m], the transpose of A
    %   b: a vector of size [m, 1]
    %   c: a vector of size [n, 1]
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
    %  rescale.left_scale_A: a vector of size  [m, 1]
    %  rescale.left_scale_F: a vector of size  [k, 1]
    %                       for rows that in the same secong order cone block, they have the same scaling factor
    %  rescale.right_scale: a vector of size [n, 1]
    %
    % explanations on some intermediate variables
    % row_norm :  store the norm of each row in  [A; F], size = [m + k, 1]
    % col_norm :  store the norm of each column in [A; F], size = [n, 1]

    t0 = tic;
    if nargin < 4
        option = 2;
    end
    % fprintf("[before rescale] A norm: %e F norm: %e \n", norm(A, 'fro'), norm(F, 'fro'));
    rescale = struct; 
    [n, m] = size(At);
    left_scale = ones(m, 1);
    right_scale = ones(n, 1);


    n_iter = 5; % number of iterations. For large model, each iteration is time consuming. Hence we only use a few iteration .



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
        elseif option == 2 % use 2 norm
            [col_norm, row_norm] = row_col_norm(At, 2, 2);
        elseif option == 3 % use infty norm
            [col_norm, row_norm] = row_col_norm(At, 3, 3);
        end
        temp_left = 1 ./ sqrt(row_norm);
        temp_right = 1 ./ sqrt(col_norm);
        At = spdiag(temp_right) * At * spdiag(temp_left);

        left_scale = left_scale .* temp_left;
        right_scale = right_scale .* temp_right;
    end

    rescale.left_scale = left_scale;
    rescale.right_scale = right_scale;

    b = b .* rescale.left_scale;
    % x = 1./ right_scale .* x;
    c = c .* rescale.right_scale;
    lb = lb ./ rescale.right_scale;
    ub = ub ./ rescale.right_scale;

    % fprintf("[after rescale]  A norm: %e F norm: %e \n", norm(A, 'fro'), norm(F, 'fro'));
    rescale.time = toc(t0);

end

