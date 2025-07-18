%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-11-27 12:28:34
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-27 20:28:38
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = Q(X, Y, cone_size)
    % if Y == 'mat', compute Q_X; else compute Q_X * Y
    detX = socp.det(X, cone_size);
    J = socp.J(cone_size);
    if strcmp(Y, 'mat')
        out = - spdiag(repelem(detX, cone_size, 1) .* J) + 2 * blk_spdiag(X, cone_size) * blk_spdiag(X, cone_size)';
    else % Y is a vector of same size as X
        % check size: to be added
        out = - repelem(detX, cone_size, 1) .* J .* Y + 2 * X .* blk_sum(X .* Y, cone_size, 1, 1) ;
    end
end