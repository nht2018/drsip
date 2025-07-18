%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-11-27 12:28:34
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-27 20:27:33
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = Qinv(X, Y, cone_size)
    % if Y == 'mat', compute Q_{X^{-1}}; else compute Q_{X^{-1}} * Y . Note that Q_{X^{-1}} = Q_{X}^{-1}
    detX = socp.det(X, cone_size);
    Xinv = socp.inv(X, cone_size);
    J = socp.J(cone_size);
    if strcmp(Y, 'mat')
        out = - spdiag(repelem(1 ./ detX, cone_size, 1) .* J) + 2 * blk_spdiag(Xinv, cone_size) * blk_spdiag(Xinv, cone_size)';
    else % Y is a vector of same size as X
        % check size: to be added
        out = - repelem(1 ./ detX, cone_size, 1) .* J .* Y + 2 * Xinv .* blk_sum(Xinv .* Y, cone_size, 1, 1) ;
    end
end