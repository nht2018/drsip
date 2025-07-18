%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-11-27 12:51:00
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-27 14:57:17
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = sqrt(X, cone_size)
    % compute X^{1/2} 
    n = size(X, 1);
    sqrtdetX = mexcone_q_ops(X, 'sqrtdet', cone_size);
    ind_head = cumsum(cone_size) - cone_size + 1;
    out = X ./ sqrt(repelem(2 * (X(ind_head) + sqrtdetX), cone_size, 1 ) );
    out(ind_head) = sqrt(0.5 * (X(ind_head) + sqrtdetX) );
end