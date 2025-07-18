%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-11-27 12:41:01
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-27 20:19:05
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = Arwinv(X, Y, cone_size)
    % Arw(X)^{-1} * Y
    % where Arw(X) = det(X) / X0 * X^{-1} * X^{-T} + 1 / X0 * diag(0, I_{n - 1})
    n = size(X, 1);
    detX = socp.det(X, cone_size);
    Xinv = socp.inv(X, cone_size);
    ind_head = cumsum(cone_size) - cone_size + 1;
    temp = ones(n, 1);
    temp(ind_head) = 0;
    if strcmp(Y, 'mat')
        out = spdiag(temp ./ repelem(X(ind_head), cone_size, 1)) + blk_spdiag(Xinv, cone_size) * spdiag(detX ./ X(ind_head)) * blk_spdiag(Xinv, cone_size)';
    else
        out = temp ./ repelem(X(ind_head), cone_size, 1) .* Y + repelem(detX ./ X(ind_head), cone_size, 1) .* Xinv .* blk_sum(Xinv .* Y, cone_size, 1, 1) ;
    end
end


% %% check 
% cone_size = [3 4 5 6]
% n = sum(cone_size);
% X = randn(n, 1);
% norm(socp.Arw(X, 'mat', cone_size) * socp.Arwinv(X, 'mat', cone_size) - speye(10), 'fro')
