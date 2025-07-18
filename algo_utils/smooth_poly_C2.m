%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-01-26 21:55:44
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-29 20:45:15
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [X, DX] = smooth_poly_C2(Z, mu, K)
% smoothing projection onto the cone K

if mu == 0
    if nargout == 1
        X = projection(Z, K);
    else
        [X, DX] = projection(Z, K);
    end
    return
end

X = MatCell(length(Z));
if nargout > 1
    DX = StructCell(length(Z));
end

for p= 1: length(K)
    cone = K{p};
    if iscell(mu)
        mup = mu{p};
    else
        mup = mu;
    end
    if strcmp(cone.type, 'l') || endsWith(cone.type, 'l')
        X{p} = mexsmooth_poly_C2(full(Z{p}), mup);
        if nargout > 1
            DX{p}.coeff = 0;
            DX{p}.shift = mexDsmooth_poly_C2(full(Z{p}), mup);
        end
    elseif strcmp(cone.type, 's')
        assert(numel(mup) == 1); % to be modified for multiple blocks
        [V, D] = eig(Z{p});  % for multiple blocks, to be optimized
        d = diag(D);
        phid = mexsmooth_poly_C2(d, mup);
        X{p} = V * spdiag(phid) * V';
        if nargout > 1
            Dphid = mexDsmooth_poly_C2(d, mup);
            DX{p}.shift = 0;
            DX{p}.coeff = mexDspectral(d, phid, Dphid);
            DX{p}.lr = V;  % save the eigenvectors
        end
    elseif strcmp(cone.type,'q') || endsWith(cone.type, 'q')
        error('not implemented yet');
        % if nargout == 1
        %     if formula == 1
        %         [X{p}] = mexprox_cone_q(full(Z{p}), mup, cone.size);  % MEX implemention
        %     else
        %         [X{p}] = mexprox_cone_q2(full(Z{p}), mup, cone.size);
        %     end
        % else
        %     if formula == 1
        %         % new version: not stable. we will dig out the reason later.
        %         [X{p}, dd, Dsch1, Dsch2, P1, P2, shift] = mexprox_cone_q(full(Z{p}), mup, cone.size);
        %         DX{p}.shift = shift;
        %         DX{p}.coeff = [Dsch1; Dsch2];
        %         DX{p}.lr = [P1, P2] ;
        %     else
        %         % old version performs more stable.
        %         % This version is compatible with Jac_lmut but not Jac_ops. Hence not work for sigma > 0 case.
        %         assert(isscalar(mup));
        %         [X{p}, DX{p}.shift, DX{p}.coeff, DX{p}.lr] = mexprox_cone_q2(full(Z{p}), mup, cone.size, det_reg);  % MEX implemention
        %     end 
        % end
    elseif strcmp(cone.type,'u')
        X{p} = full(Z{p});
        if nargout > 1
            DX{p}.coeff = 0;
            DX{p}.shift = ones(size(Z{p}));
        end
    end
end

end