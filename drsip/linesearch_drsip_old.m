%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-01-30 11:22:22
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [stepsize, X_new, Z_new, model, algo] = linesearch(X, Z, dZ, mu, model, algo, params)    
    t0 = tic;
    b = model.b;

    if ~isfield(algo, 'normA')
        algo.normA = norm(model.At);
    end
    % compute_energy = @(X_, Z_) sqrt(norm(algo.funcF(X_, Z_))^2 + mu^2);
    % compute_energy = @(X_, Z_) sqrt(norm(algo.funcF(X_, Z_))^2 + mu^2);   
    % compute_energy = @(X_, Z_) sqrt(norm(algo.AXmap(X_)-b) ^2 + mu^2 * algo.normA^2);
    compute_energy = @(X_, Z_) norm(algo.AXmap(X_)-b);
    energy_old = compute_energy(X, Z);
    energy_best = inf;
    stepsize = params.stepsize_init;
    best_stepsize = stepsize;
    while true
        Z_new = Z + dZ * stepsize;
        X_new = prox_neglog(Z_new, (1 - stepsize) * algo.AL_penalty * mu, model.K, params.det_reg);
        energy_temp = compute_energy(X_new, Z_new); 
        if energy_temp < (1 - 0.001 * stepsize) * energy_old || stepsize <= 1e-5
            best_stepsize = stepsize;
            break;
        else
            stepsize = stepsize * 0.95;
        end
    end
    stepsize = best_stepsize;


    % %% X. Chen, L. Qi, D. Sun. Global and superlinear convergence of the smoothing newton method and its application to general box constrained variational inequalities
    % funcF = @(X_, Z_) norm(algo.funcF(X_, Z_));

    % sigma = 0.01; % coeffient of the Armijo-type condition

    % F_old = norm(funcF(X, Z));
    % X0 = prox_neglog(Z, 0, model.K); % at mu=0
    % F0_old = norm(funcF(X0, Z));
    % stepsize = 0.5;
    % while true
    %     Z_new = Z + dZ * stepsize;
    %     X_new = prox_neglog(Z_new, (1 - stepsize) * algo.AL_penalty * mu, model.K);
    %     F_temp = funcF(X_new, Z_new); 
    %     if F_temp ^2 - F_old ^2 <= 2 * sigma * stepsize * F0_old^2
    %         break;
    %     elseif stepsize <= 1e-5
    %         break;
    %     else
    %         stepsize = stepsize * 0.95;
    %     end
    % end

    

    algo.linesearch.time = toc(t0);
    algo.linesearch.total_time = algo.linesearch.total_time + algo.linesearch.time;
end

