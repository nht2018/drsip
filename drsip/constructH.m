%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-01-09 20:54:38
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-04 19:20:41
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [algo] = constructH(X, y, S, mu, model, algo, params)
%% H = (I - V)^{-1} * (gamma * V + sigma * I) = - gamma * I + (sigma + gamma) *(I - V)^{-1}
K = model.K;
sigma = algo.newton.sigma; % regularization parameter for newton system
gamma = algo.AL_penalty; % augmented lagrangian penalty parameter
gammu = gamma * mu;

% if sigma > 0
if true
    invImV = algo.newton.invImV;
    % H = Jac_ops(invImV, K, 'affine', gamma, - gamma);
    H = Jac_ops(invImV, K, 'affine', sigma + gamma, - gamma);
    if strcmp(algo.newton.scaling, 'HKM')
        for p = 1:length(K)
            cone = K{p};
            if strcmp(cone.type,'q') || strcmp(cone.type, 'r2q')
                q2 = socp.inv(S{p}, cone.size);
                H{p}.lr = [socp.Qsqrt(q2, H{p}.lr(:, 1), cone.size), socp.Qsqrt(q2, H{p}.lr(:, 2), cone.size), q2];
                H{p}.coeff = [H{p}.coeff; 2 * H{p}.shift];
                H{p}.shift = H{p}.shift .* socp.det(q2, cone.size);
                H{p}.shift = repelem( - H{p}.shift, cone.size, 1) .* socp.J(cone.size);
            end
        end
    elseif strcmp(algo.newton.scaling, 'none')
        
    else
        error('unknown scaling method');
    end
else % for sigma = 0, we can compute H = gamma * (I - V)^{-1} * V = gamma * (V^{-1} - I) directly
    H = StructCell(length(K));
    proxZ = algo.newton.proxZ;
    maxH = 1e16;
    minH = 1e-16;
    for p = 1:length(K)
        cone = K{p};
        if params.vector_mu
            gammup = gammu{p};
        else
            gammup = gammu;
        end
        if strcmp(cone.type,'l') || strcmp(cone.type, 'b2l') || strcmp(cone.type, 'u2l')
            H{p}.shift = gamma * proxZ{p} .^2 ./ gammup;
            H{p}.shift = min(H{p}.shift, maxH);
            H{p}.shift = max(H{p}.shift, minH);
        elseif strcmp(cone.type, 's')
            % H.lmut{p} = @(r_) 1 / ( (1 + sigma) * mu) * proxZ{p} * r_ * proxZ{p}; % here r_ is a matrix
            error("not implemented yet");
        elseif strcmp(cone.type,'q') || strcmp(cone.type,'r2q')
            if strcmp(algo.newton.scaling, 'HKM')
                detproxZ = socp.det(proxZ{p}, cone.size);
                detproxZ = max(detproxZ, 1e-16); % avoid divide 0
                q2 = socp.inv(S{p}, cone.size);
                detq2 = soc_ops(q2, 'det', cone.size);
                H{p}.shift = gamma * (repelem( - detproxZ .* detq2 ./ gammup, cone.size, 1)) .* socp.J(cone.size) ;
                H{p}.coeff = gamma * 2 * detproxZ ./ gammup;
                H{p}.lr = socp.Qsqrt(q2, repelem (sqrt( 1 ./ detproxZ), cone.size, 1) .* proxZ{p}, cone.size);
            elseif strcmp(algo.newton.scaling, 'none')
                detproxZ = soc_ops(proxZ{p}, 'det', cone.size);
                detproxZ = max(detproxZ, 1e-16); % avoid divide 0
                H{p}.shift =  gamma * (repelem( - detproxZ ./ gammup, cone.size, 1)) .* socp.J(cone.size) ;
                H{p}.coeff = gamma * 2 * detproxZ ./ gammup;
                H{p}.lr = sqrt( 1 ./ repelem(detproxZ, cone.size, 1)) .* proxZ{p};
            else
                error('unknown scaling method');
            end
        elseif strcmp(cone.type,'u')
            H{p}.shift = 0;
        else
            error('unknown cone type');
        end
    end

end


%% handle box constraints
H_rd = H;
for p = 1:length(K)
    cone = K{p};
    if strcmp(cone.type, 'b2l')
        % reduce the dim from cone.size to cone.size / 2
        h1 = H{p}.shift(1: cone.size / 2);
        h2 = H{p}.shift(cone.size / 2 + 1: cone.size);
        H_rd{p}.shift = h1 .* h2 ./ (h1 + h2);
    end
end


algo.newton.H = H;
algo.newton.H_rd = H_rd;

maxH = 0;
minH = inf;
for p = 1:length(K)
    if ~ strcmp(K{p}.type, 'u')
        maxH = max(maxH, max(abs(H{p}.shift)));
        minH = min(minH, min(abs(H{p}.shift)));
    end
end

algo.newton.condH = maxH / minH;

end