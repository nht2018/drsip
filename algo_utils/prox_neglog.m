%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-01-26 21:55:44
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-04 19:21:57
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [X, DX] = prox_neglog(Z, mu, K, det_reg, formula)
% the proximal mapping and its derivative wrt z
% mu can be either a scalar or a cell of the same length as the cone K
% mu = 0 is equivalent to the projection mapping
% output
% X: the proximal mapping of Z
% DX: the derivative of the proximal mapping wrt Z
%     DX.lmut is a functional handle, i.e., DX.lmut(Z_) = DX.lmut * Z_

if nargin < 4
    det_reg = 0;
end
if nargin < 5
    formula = 2;
end
assert(ismember(formula, [1, 2]));

if mu == 0
    if nargout == 1
        X = projection(Z, K);
    else
        [X, DX] = projection(Z, K);
    end
    return
end

X = MatCell(length(Z));
phi =  @(z, mu) 0.5 * (z + sqrt(z.^2 + 4 * mu));
if nargout > 1
    DX = StructCell(length(Z));
    Dphi = @(z, mu) 0.5 * ( 1 + z ./ sqrt(z .^ 2 + 4 * mu) );
end

for p= 1: length(K)
    cone = K{p};
    if iscell(mu)
        mup = mu{p};
    else
        mup = mu;
    end
    if strcmp(cone.type, 'l') || endsWith(cone.type, 'l')
        X{p} = phi(full(Z{p}), mup);
        if nargout > 1
            DX{p}.coeff = 0;
            DX{p}.shift = Dphi(full(Z{p}), mup);
        end
    elseif strcmp(cone.type, 's')
        assert(numel(mup) == 1); % to be modified for multiple blocks
        [V, D] = eig(Z{p});  % for multiple blocks, to be optimized
        d = diag(D);
        phid = phi(d, mup);
        X{p} = V * spdiag(phid) * V';
        if nargout > 1
            Dphid = Dphi(d, mup);
            DX{p}.shift = 0;
            DX{p}.coeff = mexDspectral(d, phid, Dphid);
            DX{p}.eigs = V;  % save the eigenvectors
        end
    elseif strcmp(cone.type,'q') || endsWith(cone.type, 'q')
        if nargout == 1
            if formula == 1
                [X{p}] = mexprox_cone_q(full(Z{p}), mup, cone.size);  % MEX implemention
            else
                [X{p}] = mexprox_cone_q2(full(Z{p}), mup, cone.size);
            end
        else
            if formula == 1
                % new version: not stable. we will dig out the reason later.
                [X{p}, dd, Dsch1, Dsch2, P1, P2, shift] = mexprox_cone_q(full(Z{p}), mup, cone.size);
                DX{p}.shift = shift;
                DX{p}.coeff = [Dsch1; Dsch2];
                DX{p}.lr = [P1, P2] ;
            else
                % old version performs more stable.
                % This version is compatible with Jac_lmut but not Jac_ops. Hence not work for sigma > 0 case.
                assert(isscalar(mup));
                [X{p}, DX{p}.shift, DX{p}.coeff, DX{p}.lr] = mexprox_cone_q2(full(Z{p}), mup, cone.size, det_reg);  % MEX implemention
            end
            
        end
    elseif strcmp(cone.type,'u')
        X{p} = full(Z{p});
        if nargout > 1
            DX{p}.coeff = 0;
            DX{p}.shift = ones(size(Z{p}));
        end
    end
end


% %% check prox_neglog by the identity Z = X - mu * inv(X)
% fprintf("check prox_neglog: %.6e\n", norm(Z - X + mu * Dneglog(X, K)) / norm(Z));

% if nargout > 1
%     % check DX.lmut by the identity DX.lmut * (2 * X - Z) == X
%     fprintf("check DX relresd = %e\n", norm(Jac_lmut(K, DX, 2 * X - Z) - X) / norm(X));
% end
end