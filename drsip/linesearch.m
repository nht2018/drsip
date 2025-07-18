%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-11 16:18:48
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [stepsize, X_new, y_new, S_new, model, algo] = linesearch(X, y, S, dX, dy, dS, mu, model, algo, params)    
    t0 = tic;
    
    energy_func = @(X_, y_, S_) norm(algo.presmap(X_))^2 + norm(algo.dresmap(y_, S_))^2 + norm(X_ - algo.smooth_func(X_ - algo.AL_penalty * S_, algo.mu_decay * algo.AL_penalty * mu))^2 ;
    % energy_func = @(X_, y_, S_) norm(algo.presmap(X_))^2 + norm(algo.dresmap(y_, S_))^2 + norm(X_ - algo.projection_func(X_ - algo.AL_penalty * S_))^2 ;

    energy_old = energy_func(X, y, S);
    energy_best = inf;
    stepsize = params.stepsize_init;
    best_stepsize = stepsize;
    while true
        X_new = X + stepsize * dX;
        y_new = y + stepsize * dy;
        S_new = S + stepsize * dS;
        energy_temp = energy_func(X_new, y_new, S_new);
        if energy_temp < energy_best
            energy_best = energy_temp;
            best_stepsize = stepsize; 
        end
        if energy_best <= energy_old || stepsize <= 1e-5
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

