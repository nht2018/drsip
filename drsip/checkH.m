%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-11-24 12:34:06
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-28 12:15:05
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   

% p = 2;
% matH = spdiag(H.shift{p}) + H.eigvec{p} * spdiag(H.coeff{p}) * H.eigvec{p}';

% sigma = 0;
% Htemp = Jac_ops(DX, K, 'affine_inv', 1, -(1+sigma));
% matHtemp = spdiag(repelem(Htemp.shift{p}, K{p}.size, 1)) + Htemp.eigvec{p} * spdiag(Htemp.coeff{p}) * Htemp.eigvec{p}';
% H1 = Jac_ops(Htemp, K, 'affine', -(1+2*sigma), -1);
% matH1 = spdiag(repelem(H1.shift{p}, K{p}.size, 1)) + H1.eigvec{p} * spdiag(H1.coeff{p}) * H1.eigvec{p}';


% matV = spdiag(repelem(DX.shift{p}, K{p}.size, 1)) + DX.eigvec{p} * spdiag(DX.coeff{p}) * DX.eigvec{p}';

% check H^{-1} = V^{-1} - I by checking (I - V) * H1 = V
% norm((speye(size(matV)) - matV) * matH - matV, 'fro')



% % check Htemp = (DX-I)^{-1}
% norm((matV - speye(size(matV)) ) * matHtemp - speye(size(matV)), 'fro')


% norm(matH - matH1, 'fro')


% % check Htemp = (DX-I)^{-1}
% r = MatCell.rand_like(Z);
% r1 = Jac_lmut(K, Htemp, r);
% r2 = Jac_lmut(K, DX, r1) - r1;
% norm(r2 - r)

% % check H1^{-1} = DX^{-1} - I by checking (I - DX) * H1 = DX
% r = MatCell.rand_like(Z);
% r1 = Jac_lmut(K, H1, r);
% r2 = r1 - Jac_lmut(K, DX, r1);
% r3 = Jac_lmut(K, DX, r);
% norm(r2 - r3)

% check H^{-1} = V^{-1} + I by checking (I + V) * H = V
% r = MatCell.rand_like(Z);
% r1 = Jac_lmut(K, H, r);
% r2 = r1 + Jac_lmut(K, V, r1);
% r3 = Jac_lmut(K, V, r);
% norm(r2 - r3)

% % check matH = 1 / mu (- detx * J + 2 X X^T)
% cone = K{p};
% ind_head = cumsum(cone.size) - cone.size + 1;  % record the start index of each cone blocks
% tempdiag = ones(size(Z{p}));
% tempdiag(ind_head) = -1;
% res = matH1 - (1/mu) * ( spdiag(repelem( detx, cone.size, 1) .* tempdiag) + 2 * blk_spdiag(X{p}, cone.size) * blk_spdiag(X{p}, cone.size)');

% % norm(spdiag(H.shift{p}) - (1/mu) * spdiag(repelem( detx, cone.size, 1) .* tempdiag), 'fro')
% % norm( H.eigvec{p} * spdiag(H.coeff{p}) * H.eigvec{p}' - (1/mu) * 2 * blk_spdiag(X{p}, cone.size) * blk_spdiag(X{p}, cone.size)', 'fro' )


% %check DX^{-1} - I = mu * ( 1 / detx * J + 2 Xinv Xinv^T)
% %  by checking I - DX = mu * ( 1 / detx * J + 2 Xinv Xinv^T) * DX
% Xinv = soc_ops(X{p}, 'inv', cone.size);
% % check X - mu * Xinv = Z
% norm(X{p} - mu * Xinv - Z{p}, 'fro')
% % norm(speye(size(matV)) - matV - ...
% % mu * ( spdiag(repelem( 1 ./ detx, cone.size, 1) .* tempdiag) + 2 * blk_spdiag(Xinv, cone.size) * blk_spdiag(Xinv, cone.size)') * matV, 'fro')


% check lhs = A * H * A'
y0 = randn(size(y));
l1 = algo.AXmap(Jac_lmut(K, H, algo.Atymap(y0)));
l2 = lhs.lmut(y0);
norm(l1 - l2)
