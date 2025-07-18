%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-01-09 20:54:38
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-10 16:18:18
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [algo] = constructH(X, y, S, mu, model, algo, params)
%% H = ((1+tau_xs)I - V)^{-1} * (gamma * V) = gamma * (1+tau_xs) * ((1+tau_xs)I - V)^{-1} - gamma * I
K = model.K;
gamma = algo.AL_penalty; % augmented lagrangian penalty parameter
gammu = gamma * mu;

tau_xs = algo.newton.tau_xs; % regularization parameter for newton system
invImV = algo.newton.invImV;
% H = Jac_ops(invImV, K, 'affine', gamma, - gamma);
H = Jac_ops(invImV, K, 'affine',  gamma * (1+tau_xs), - gamma);
if strcmp(algo.newton.scaling, 'HKM')
    for p = 1:length(K)
        cone = K{p};
        if strcmp(cone.type,'q') || endsWith(cone.type, 'r2q')
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