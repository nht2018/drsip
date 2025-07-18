%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-01-26 21:55:44
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-11 12:23:03
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
%%*****************************************************************
%% drs_dir: solve the Newton equation and obtain the direction
%%*****************************************************************

function [dX, dy, dS, res_rd, res_nt, model, algo] = sn_dir(X, y, S, mu, model, algo, params, pred_or_corr, dX_pred, dS_pred)

assert(ismember(pred_or_corr, {'pred', 'corr'}), 'pred_or_corr should be either pred or corr');
At = model.At;
K = model.K;

gamma = algo.AL_penalty; % augmented lagrangian penalty parameter
gammu = gamma * mu;
% gammu = mu;


if ~isfield(model, 'idx_u')
    model.idx_u = []; % block index of unbounded cone
    model.size_u = [];
    for p = 1: length(K)
        cone = K{p};
        if strcmp(cone.type, 'u')
            model.idx_u = [model.idx_u, p];
            model.size_u = [model.size_u, cone.size];
        end
    end
end
algo.newton.q_dense_col = 0;

%% scaling the direction
if strcmp(pred_or_corr, 'pred')
    tildeX = X;
    utildeS = S;
    if strcmp(algo.newton.scaling, 'HKM')
        for p = 1: length(K)
            cone = K{p};
            if strcmp(cone.type, 'q')
                p2 = S{p};
                tildeX{p} = socp.Qsqrt(p2, X{p}, cone.size);
    %             utildeS{p} = socp.Qsqrtinv(p2, S{p}, cone.size);
                utildeS{p} = socp.e(cone.size);
                assert(norm(utildeS{p} - socp.e(cone.size)) < 1e-10, 'utildeS is not unit');
            end
        end
    elseif strcmp(algo.newton.scaling, 'none')

    else
        error('unknown scaling method');
    end
    Z = tildeX - gamma * utildeS;
    [Phiz, V] = algo.smooth_func(Z, gammu);

    algo.newton.Z = Z;
    algo.newton.Phiz = Phiz;
    algo.newton.tildeX = tildeX;
    algo.newton.utildeS = utildeS;

    %% compute rhs in the reduced system
% dX + H * dS = r_c = - (I-V)^{-1} (X - Phiz)             complementarity
% A * dX = r_p = - (A x - b )                   primal feasibility
% A' * dy + dS = r_d = - (A' y + S - c)         dual feasibility
    r_p = - algo.presmap(X) ;
    r_d = - algo.dresmap(y, S);
    [Phiz1] = algo.smooth_func(Z, algo.mu_decay * gammu);
    r_xs = Phiz1 - tildeX;

    %% compute direction
    if params.newton_reg
        algo.newton.tau_xs = max(algo.newton.kappa * (norm(r_xs) ^ algo.newton.eta), 1e-12);
        algo.newton.tau_p = algo.newton.kappa * (norm(r_p) ^ algo.newton.eta);
        algo.newton.tau_d = algo.newton.kappa * (norm(r_d) ^ algo.newton.eta);
    else
        algo.newton.tau_xs = 0;
        algo.newton.tau_p = 0;
        algo.newton.tau_d = 0;
    end


    invImV = Jac_ops(V, K, 'affine_inv', - 1, 1 + algo.newton.tau_xs); % ((1+ tau_xs)I - V)^{-1}
    algo.newton.V = V;
    algo.newton.invImV = invImV;
    r_c = Jac_lmut(K, invImV, r_xs);

    if strcmp(algo.newton.scaling, 'HKM')
        for p = 1: length(K)
            cone = K{p};
            if strcmp(cone.type, 'q')
                p2 = S{p};
                r_c{p} = socp.Qsqrtinv(p2, r_c{p}, cone.size);
            end
        end
    end

    if ~isempty(model.idx_u)
        r_c{model.idx_u} = zeros(sum(model.size_u), 1);
    end

    algo.newton.r_p = r_p;
    algo.newton.r_d = r_d;
    algo.newton.r_xs = r_xs;
    algo.newton.r_c = r_c;
    t = tic;
    algo = constructH_SN(X, y, S, mu, model, algo, params);
    H = algo.newton.H;
    H_rd = algo.newton.H_rd;

    [lhs, model, algo] = constructlhs(H_rd, model, algo, params);
    algo.newton.reg_init = (1+algo.newton.tau_d) * algo.newton.tau_p;
    % if algo.iter == 1
    %     fprintf(algo.fid, "[newton system] dim1 = %d, dim2 = %d, dim3 = %d\n", lhs.dim1, lhs.dim2, lhs.dim3);
    %     if algo.newton.system_opt == 1 || algo.newton.system_opt == 3
    %         fprintf(algo.fid, "\tnonzeros = %d, sparsity = %f\n", nnz(lhs.mat), nnz(lhs.mat) / numel(lhs.mat));
    %     end
    % end
    algo.newton.lhs = lhs;
    algo.newton.time_construct = toc(t);
else
    tildeX = algo.newton.tildeX;
    utildeS = algo.newton.utildeS;
    Z = algo.newton.Z;
    Phiz = algo.newton.Phiz;
    V = algo.newton.V;
    invImV = algo.newton.invImV;

    r_p = algo.newton.r_p;
    r_d = algo.newton.r_d;
    % [Phiz1] = algo.smooth_func(Z, algo.mu_decay * gammu);
    % tilde_dX_pred = dX_pred;
    % utilde_dS_pred = dS_pred;
    % [Phiz2] = algo.smooth_func(Z + tilde_dX_pred - gamma * utilde_dS_pred, algo.mu_decay * gammu);
    % r_xs = Phiz1 - tildeX + Phiz2 - tildeX - tilde_dX_pred;
    % r_c = Jac_lmut(K, invImV, r_xs);

    %  linearize proximal
    [Phiz1] = algo.smooth_func(Z, algo.mu_decay * gammu);
    r_xs = Phiz1 - tildeX;
    for p =1: length(K)
        cone = K{p};
        if strcmp(cone.type, 'l') || strcmp(cone.type, 'b2l')
            temp = max(abs(dX_pred{p} - dS_pred{p}), sqrt(Z{p} .^ 2 + 4 * algo.mu_decay * gammu) ); 
            temp = max(temp, 1);
            corr = 0.5 * (dX_pred{p} -  dS_pred{p}) .^ 2 .* (2 * algo.mu_decay * gammu) ./ ( temp .^(3) );
            corr = corr * 0.1;
            ratio = 1;
            bound = ratio * r_xs{p};
            idx = abs(corr) >= bound; 
            corr(idx) = sign(corr(idx)) .* bound(idx);
            % fprintf("l: corr / r_xs = %f\n", max(abs(corr ./ r_xs{p})))
            r_xs{p} = r_xs{p} + corr;
        elseif strcmp(cone.type, 'q')
            ind_head = cumsum(cone.size) - cone.size + 1;
            eigz1 = Z{p}(ind_head) + socp.sqrtdet(Z{p}, cone.size);
            eigz1 = max(eigz1, 1);
            eigz2 = Z{p}(ind_head) - socp.sqrtdet(Z{p}, cone.size);
            eigz2 = max(eigz2, 1);

            u = dX_pred{p} - algo.AL_penalty * dS_pred{p};
            v1 = Z{p} ./ repelem(max(socp.norm_xbar(Z{p}, cone.size), 1e-16) , cone.size, 1);
            v1(ind_head) = 1;
            v2 = - v1;
            v2(ind_head) = 1;
            corr = repelem((2 * algo.mu_decay * gammu) ./ sqrt(eigz1 .^ 2 + 4 * algo.mu_decay * gammu) .* socp.inner_product(u, v1, cone.size) .^ 2, cone.size, 1) .* v1 + ...
                repelem((2 * algo.mu_decay * gammu) ./ sqrt(eigz2 .^ 2 + 4 * algo.mu_decay * gammu) .* socp.inner_product(u, v2, cone.size) .^ 2, cone.size, 1) .* v2;

            corr = corr * 0.1;
            ratio = 1;
            bound = ratio * r_xs{p};
            idx = abs(corr) >= bound; 
            corr(idx) = sign(corr(idx)) .* bound(idx);
            % fprintf("q: corr / r_xs = %f\n", max(abs(corr ./ r_xs{p})))
            r_xs{p} = r_xs{p} + corr;
        end
    end
    r_c = Jac_lmut(K, invImV, r_xs);

    if strcmp(algo.newton.scaling, 'HKM')
        for p = 1: length(K)
            cone = K{p};
            if strcmp(cone.type, 'q')
                p2 = S{p};
                r_c{p} = socp.Qsqrtinv(p2, r_c{p}, cone.size);
            end
        end
    end

    if ~isempty(model.idx_u)
        r_c{model.idx_u} = zeros(sum(model.size_u), 1);
    end

    
    H = algo.newton.H;
    H_rd = algo.newton.H_rd;
    lhs = algo.newton.lhs;
end





%% handle box constraints
r_c_rd = r_c;
dS_shift = MatCell(length(K));
for p = 1:length(K)
    cone = K{p};
    if strcmp(cone.type, 'b2l')
        % reduce the dim from cone.size to cone.size / 2
        h1 = H{p}.shift(1: cone.size / 2);
        h2 = H{p}.shift(cone.size / 2 + 1: cone.size);
        rbox = X{p}(1: cone.size / 2) + X{p}(cone.size / 2 + 1: cone.size) - cone.params(:, 2) + cone.params(:, 1);
        dS_shift{p} = (r_c{p}(1: cone.size / 2) + r_c{p}(cone.size / 2 + 1: cone.size) - rbox) ./ (h1 + h2);
        r_c_rd{p} = r_c{p}(1: cone.size / 2) - h1 .* dS_shift{p};
    end
end

% rhs_rd = r_p + A * ( H * r_d - r_c)
rhs_rd = (1+algo.newton.tau_d) * r_p + AXfun(model.K, model.At, (Jac_lmut(K, H_rd, r_d) - (1+algo.newton.tau_d) * r_c_rd));
if ~isempty(model.idx_u)
    rhs_u = r_d{model.idx_u};
else
    rhs_u = [];
end

if algo.newton.system_opt == 1
    rhs = [rhs_rd;
        rhs_u;
        zeros(lhs.dim3, 1)];
else % algo.newton.system_opt == 3 || algo.newton.system_opt == 4
    rhs = [rhs_rd;
        rhs_u];
end
rhs = full(rhs) ;


%% factorize the reduced system
if strcmp(pred_or_corr, 'pred')
    t = tic;
    if strcmp(algo.newton.system_opt, 1) % to be improved
        algo.newton.ordering = [];
    end
    [algo.newton] = linsysfact(lhs, algo.newton);
    algo.newton.time_fact = toc(t);
    algo.newton.total_time_fact = algo.newton.total_time_fact + algo.newton.time_fact;
end


%% solve the reduced system
t = tic;
[w, algo.newton] = linsyssolve(lhs, rhs, algo.newton);
% w = lhs.mat  \ rhs;

%% check residual of reduced system
if  any(isnan(w)) 
    error("w is nan\n");
end
if  any(isinf(w)) 
    error("w is inf\n");
end
res_rd = norm(rhs - lhs.lmut(w)) / (1 + norm(rhs));

%% recover variables in Newton system
%% compute w
dy = w(1: lhs.dim1); % if preconditioner is used, then d(1: lhs.dim1) actually is R^{-1} * w
dX_u = w(lhs.dim1 + 1: lhs.dim1 + lhs.dim2); 

%% compute dX and dS
dS = r_d - algo.Atymap(dy);
%% handle box constraints
for p = 1:length(K)
    cone = K{p};
    if strcmp(cone.type, 'b2l')
        % reduce the dim from cone.size to cone.size / 2
        h1 = H{p}.shift(1: cone.size / 2);
        h2 = H{p}.shift(cone.size / 2 + 1: cone.size);
        dS{p} = [h2 ./ (h1 + h2) .* dS{p} + dS_shift{p};
            - h1 ./ (h1 + h2) .* dS{p} + dS_shift{p}];
    end
end

dX = r_c - Jac_lmut(K, H, dS);


if ~isempty(model.idx_u)
    dX{model.idx_u} = dX_u;
    dS{model.idx_u} = zeros(sum(model.size_u), 1);
end

for p = 1: length(K)
    cone = K{p};
    if ~ strcmp(cone.type, 's')
        if issparse(dX{p})
            dX{p} = full(dX{p});
        end
        if issparse(dS{p})
            dS{p} = full(dS{p});
        end

    end
end

% %% check box constrains
% for p = 1: length(K)
%     cone = K{p};
%     if strcmp(cone.type, 'b2l')
%         fprintf("check box: %.4e\n", norm(dX{p}(1: cone.size / 2) + dX{p}(cone.size / 2 + 1: end)) );
%     end
% end

%% check residual of Newton system
if params.print_log == 2
    tau_p = algo.newton.tau_p;
    tau_d = algo.newton.tau_d;
    tau_xs = algo.newton.tau_xs;
    res_p = norm(algo.AXmap(dX) + tau_p * dy - r_p) ;
    res_d = norm(algo.dresmap(dy, (1+tau_d) * dS) + model.c - r_d);

    % res_c = (I - V) * dX + gamma * V * dS + X - Phiz1
    if strcmp(algo.newton.scaling, 'HKM')
        tilde_dX = dX;
        utilde_dS = dS;
        for p = 1: length(K)
            cone = K{p};
            if strcmp(cone.type, 'q')
                p2 = S{p};
                tilde_dX{p} = socp.Qsqrt(p2, tilde_dX{p}, cone.size);
                utilde_dS{p} = socp.Qsqrtinv(p2, utilde_dS{p}, cone.size);
            end
        end
        res_c = (1+tau_xs)*tilde_dX - Jac_lmut(K, V, tilde_dX) + gamma * Jac_lmut(K, V, utilde_dS) - r_xs;
    elseif strcmp(algo.newton.scaling, 'none')
        res_c =  (1+tau_xs)*dX - Jac_lmut(K, V, dX) + gamma * Jac_lmut(K, V, dS) - r_xs;
    end

    % res_c = dX + Jac_lmut(K, H, dS) - r_c;
    if ~isempty(model.idx_u)
        res_c{model.idx_u} = zeros(sum(model.size_u), 1);
    end
    res_c = norm(res_c);
    res_nt = sqrt(norm(res_p) ^2 + norm(res_d) ^2 + norm(res_c) ^2) / (1 + sqrt(norm(r_p) ^2 + norm(r_d) ^2 + norm(r_c) ^2));
end


algo.newton.time_solve = toc(t);
algo.newton.total_time_solve = algo.newton.total_time_solve + algo.newton.time_solve;


end