%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-23 17:20:41
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-23 18:00:17
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   

function out = Qsqrtinv(X, Y, cone_size)
    % if Y == 'mat', compute Q_{X^{1/2}}; else compute Q_{X^{1/2}} * Y
    if strcmp(Y, 'mat')
        % Q_{X^{1/2}} = [X0 , Xbar;
        %             Xbar,  sqrt(detX) * I + 1/(X0 + sqrt(detX)) * Xbar * Xbar']
        error('not implemented yet')
    else % Y is a vector of same size as X
        ind_head = cumsum(cone_size) - cone_size + 1;
        sqrtdetX = socp.sqrtdet(X, cone_size);
        Jinprod = socp.J_inner_product(X, Y, cone_size);
        temp = - (Y(ind_head) .* sqrtdetX + Jinprod) ./ (sqrtdetX + X(ind_head)); 
        out  = repelem(sqrtdetX, cone_size, 1) .* Y + repelem(temp, cone_size, 1) .* X; 
        out(ind_head) = Jinprod; 
        out = out ./ repelem(sqrtdetX .^ 2, cone_size, 1);
    end
    
end