%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-21 17:12:17
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-01-02 21:46:15
%
%  Copyright (c) 2023, Hantao Nie, Peking University.
%

% solve standard form conic programming
% min <c, x>
% st  A * x == b
%     x in K
% input model has fields:
% At: MatCell
% K: a cell of BasicCones
% c: MatCell
% b: vector
% or 
% solve general form conic programming
% min <c, x> + c0
% st  lc <= A * x == uc
%     F * x + g in D
%     lx <= x <= ux



function [out, hist] = drsip(model_in, params, init_point)

addpath(genpath("../algo_utils"));
t0 = tic;

if nargin < 3
    init_point = [];
end

%% set parameters
if ~exist('params', 'var'); params = struct; end
params = set_default_params(params);


%% check input model
if check_std_model(model_in)
    std_model = model_in;
    model = std_model;
elseif check_gen_model(model_in)
    std_model = standardize_forward(model_in);
    assert(check_std_model(std_model));
    model = std_model;
else
    error("Input model is not a standard form or general form socp");
end


model.At = reshape(model.At, [], 1);
model.c = reshape(model.c, [], 1);


% %% check linear independence and remove redundant constraints
% [model.At, model.b, ~, indeprows, depconstr, feasible, model.AAt] = ...
% checkdepconstr(Cone.toblk(model.K), model.At, model.b, zeros(size(model.b)), 1);
% if size(model.At{1}, 2) < size(std_model.At{1}, 2)
%     fprintf("Removed %d redundant constraints\n", size(std_model.At{1}, 2) - size(model.At{1}, 2));
% end


algo = struct; % record intermediate information during the algorithm


if params.warning_on; warning('on', 'all'); else; warning('off', 'all'); end

if strcmp(params.log_path, '')
    fid = 1; % fid=1: print to the console
else
    [log_dir, ~, ~] = fileparts(params.log_path);
    if ~exist(log_dir, 'dir')
        mkdir(log_dir);
    end
    fid = fopen(params.log_path, 'w');
end

algo.fid = fid;
fprintf(fid, "[%s]\n", char(datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss')));  % Display current date and time

%% record parameters
if params.print_params
    fprintf(fid, "[params]\n");
    fn = fieldnames(params);
    for i = 1:numel(fn)
        fprintf(fid, '\t%s: %s\n', fn{i}, num2str(params.(fn{i})));
    end
end

%% print problem informations
fprintf(fid, "Input problem\n");
print_std_model(model, fid);

%% convert rotated quadratic cone to quadratic cone
model = cone_r2q(model);


%% rescale model
if params.rescale_option ~= 0
    [model, algo] = row_col_rescale(model, algo, params.rescale_option);
    fprintf(fid, "[rescale A by row and column]\n");
    fprintf(fid, "\tmax(left scale): %e , min(left scale): %e\n", max(algo.rescale.left_scale), min(algo.rescale.left_scale));
    fprintf(fid, "\tmax(right scale): %e , min(right scale): %e\n", MatCell.max_all(algo.rescale.right_scale), MatCell.min_all(algo.rescale.right_scale));
    fprintf(fid, "\tRescaling model costs %.2f seconds\n", algo.rescale.time);
    % fprintf(fid, "Rescaled problem\n");
    % print_std_model(model, fid);
end




if params.u2l
    % % detect unbouned block hiden in linear block
    % [model, ublkinfo] =  cone_l2u(model);
    
    % convert unbounded block to linear block
    [model] = cone_u2l(model);
end

%% convert box block to linear cone
[model] = cone_b2l(model) ;

%% dense detect
if MatCell.sparsity(model.At) >= params.densemat_threshold
    fprintf(fid, "[dense detect] A is dense\n");
    model.At = MatCell.full(model.At);
    model.is_dense = true;
else
    fprintf(fid, "[dense detect] A is sparse\n");
    model.is_dense = false;
end




fprintf(fid, "Preprocessed problem\n");
print_std_model(model, fid);



%% dense column detect
algo.model.col_den = MatCell(length(model.At));
algo.model.col_sp = MatCell(length(model.At));
algo.model.At_col_den = MatCell(length(model.At));
algo.model.At_col_sp = MatCell(length(model.At));
algo.model.AHlrspcol  = MatCell(length(model.K));
algo.model.AHlrdencol = MatCell(length(model.K));
fprintf(fid, "[dense column detect]\n");
for p = 1: length(model.At)
    [algo.model.col_den{p}, algo.model.col_sp{p}] = detect_dense_row(model.At{p});
    algo.model.At_col_den{p} = model.At{p}(algo.model.col_den{p}, :);
    algo.model.At_col_sp{p} = model.At{p}(algo.model.col_sp{p}, :);
    fprintf(fid, "\tblock %d: %d dense columns, %d sparse columns\n", p, length(algo.model.col_den{p}), length(algo.model.col_sp{p}));
end



need_factAAt = true;

if need_factAAt    
    %% factorize AAt
    [model, algo] = factAAt2(model, algo, params);
    fprintf(fid, "[factorize A * A']\n");

    if isfield(algo.model, 'detect_very_sparse_rows') && algo.model.detect_very_sparse_rows
        fprintf(fid, "\tDetect %i very sparse rows in A * A'\n", numel(algo.model.row_sp));
    end
    fprintf(fid, "\tActual dimension: %d\n", algo.AAt.dim);
    if strcmp(algo.AAt.solver, 'ldlchol')
        fprintf(fid, "\tldlchol: LD factor nonzeros: %d\n", nnz(algo.AAt.LD));
    elseif strcmp(algo.AAt.solver, 'lchol')
        fprintf(fid, "\tlchol: L factor nonzeros: %d\n", nnz(algo.AAt.L));
    elseif strcmp(algo.AAt.solver, 'spchol')
        fprintf(fid, "\tspchol: R factor nonzeros: %d\n", nnz(algo.AAt.R));
    elseif strcmp(algo.AAt.solver, 'chol')
        fprintf(fid, "\tchol: R factor nonzeros: %d\n", nnz(algo.AAt.R));
    else
        fprintf(fid, "\tunknown solver: %s\n", algo.AAt.solver);
    end
    if isfield(algo.AAt, 'reg') && algo.AAt.reg > 0
        fprintf(fid, "\tAdd regularization %.2e to A * A'\n", algo.AAt.reg);
    end
    fprintf(fid, "\tUse %s to handle dense columns\n", algo.AAt.dense_column_strategy);
    fprintf(fid, "\tFactorizing  A * A' costs %.2f seconds\n", algo.AAt.time_fact);
    
end
%% precondition on A
if params.precond_A
    if model.n_box %% to be implemented in the future
        error("Preconditioning on A is not supported when there are box constraints");
    end
    %% set parameters for factorizing  A * A'
    precond = algo.AAt;
    assert(strcmp(precond.dense_column_strategy, 'none') );
    precond = rmfield(precond, 'time_fact');
    model.b = algo.AAt.fwsolve(model.b);
    model.precond = precond;
else
    model.precond = struct;
    model.precond.isidentity = true;
end



%% scale b and c
model_unscaled = model;
if params.scale_bc_flag
    if params.scale_b ~= 0
        scale_b = params.scale_b;
    else
        % scale_b = norm(model.b);
        scale_b = max(norm(model.b), 1);
    end
    if params.scale_c ~= 0
        scale_c = params.scale_c;
    else
        scale_c = max(norm(model.c), 1);
    end
    model.b = model.b / scale_b;
    model.c = model.c / scale_c;
    if model.n_box
        model.b_int = model.b_int / scale_b;
        model.c_int = model.c_int / scale_c;
        for p = 1: length(model.K)
            if strcmp(model.K{p}.type, 'b2l')
                model.K{p}.params = model.K{p}.params / scale_b;
            end
        end
    else
        model.b_int = model.b;
        model.c_int = model.c;
    end
    fprintf(fid, "[scale]\n");
    fprintf(fid, "\tscale primal variables by %e\n", scale_b);
    fprintf(fid, "\tscale dual variabls by %e\n", scale_c);
else
    scale_b = 1;
    scale_c = 1;
end
algo.scale.scale_b = scale_b;
algo.scale.scale_c = scale_c;



algo.AXmap = @(x) AXmap(x, model.K, model.At, model.precond);
% algo.AXmap_int = @(x) AXfun(model.K, model.At_int, x);
algo.AXmap_int = @(x) AXmap_int(x, model.K, model.At, model.precond);
algo.Atymap = @(y) Atymap(y, model.K, model.At, model.precond);
algo.Atymap_int = @(y) Atymap(y, model.K, model.At_int, model.precond);
% algo.Atymap_int = @(y, n_box) Atymap_int(y, model.K, model.At, model.precond, n_box);
algo.presmap = @(x) presmap(model.K, model.At, x, model.b, model.precond);
% algo.presmap_int = @(x) algo.AXmap_int(x) - model.b_int;
algo.presmap_int = @(x) presmap_int(model.K, model.At, x, model.b, model.precond);
algo.dresmap = @(y, S) dresmap(model.K, model.At, y, S, model.c, model.precond);
algo.dresmap_int = @(y, S) algo.Atymap_int(y) + S - model.c_int;
algo.Adagger = @(y) algo.Atymap_int((invAAt(y, model, algo, params)) );    % if preconditioning is used, Adagger is equivalent to Atymap


%% set parameters for newton solver and initialze
algo.newton = struct;
algo.newton.initialized = false;  % a flag to indicate whether ordering and symbolic factorization has been done, only solver 'ldlchol', 'pardiso' need this flag
res_rd_pred = 0;
res_nt_pred = 0;
res_rd_corr = 0;
res_nt_corr = 0;
algo.newton.time_construct = 0;
algo.newton.time_fact = 0;
algo.newton.time_solve = 0;
algo.newton.total_time_fact = 0;
algo.newton.total_time_solve = 0;
algo.newton.dim = 0;
algo.newton.pcg_iter = 0;
algo.newton.newton_tol = params.newton_tol;
algo.newton.q_dencol = 0;
algo.newton.sigma = 0;
algo.newton.kappa = 1e-6;
algo.newton.sigPow = 1;
algo.newton.reg = 0;
algo.newton.condH = 0;
algo.newton.method = "undefined";

%% select solver for newton system
% the automatic scheme is to be improved in the future

% check whether there exist free variables
exist_u = 0;
exist_q = 0;
exist_s = 0;
for p=1:length(model.K)
    cone = model.K{p};
    if strcmp(cone.type, 'u')
        exist_u = 1;
    end
    if strcmp(cone.type, 'q')
        exist_q = 1;
    end
    if strcmp(cone.type, 's')
        exist_s = 1;
    end
end

algo.model.exist_u = exist_u;
algo.model.exist_q = exist_q;
algo.model.exist_s = exist_s;
%% I think first decide system_opt and then decide solver is better.
%% decide strategy of constructing newton system
if params.system_opt == 0 % automatic decide
    if strcmp(params.newton_solver, 'pcg') || exist_s % pcg or semidefinite cone
        algo.newton.system_opt = 4; % iterative system to match pcg
    elseif model.is_dense || (~ exist_q && ~ exist_u) % dense system or only linear cone
        algo.newton.system_opt = 3; % dense system
    else
        algo.newton.system_opt = 1; % sparse system
    end
else % user specified
    algo.newton.system_opt = params.system_opt;
end

%% decide solver for newton system
if strcmp(params.newton_solver, 'default') % automatic decide solver
    if algo.newton.system_opt == 4 % iterative system
        algo.newton.solver = 'pcg';
    else
        algo.newton.solver = 'auto';
    end
else % user specified solver
    algo.newton.solver = params.newton_solver;
end
fprintf(fid, "[Newton system]\n");
fprintf(fid, "\tSolver for newton system: %s\n", algo.newton.solver);

assert(ismember(algo.newton.system_opt, [1, 2, 3, 4])); % now only support sparse system, dense system and iterative system
if exist_s
    assert(algo.newton.system_opt == 4 && strcmp(algo.newton.solver, 'pcg'), "Semidefinite cone can only be used with iterative system");
end
fprintf(fid, "\tReduced system option: %d\n", algo.newton.system_opt);
if algo.newton.system_opt == 1
    if isfield(algo.model, 'n_dense_col') && algo.model.n_dense_col > 0
        fprintf(fid, "\tUse %s to handle dense columns\n", params.dense_column_strategy);
    end
end
algo.newton.scaling = params.scaling_method;
fprintf(fid, "\tScaling method: %s\n", algo.newton.scaling);

% add support for pardiso
if strcmp(algo.newton.solver, 'pardiso')
    addpath("pardiso-matlab-recipes");
end


% check match between solver and reduced system option
if strcmp(algo.newton.solver, 'chol') || strcmp(algo.newton.solver, 'lchol')
    assert(algo.newton.system_opt == 3 || algo.newton.system_opt == 2, "Cholesky decomposition can only be used for dense system");
end
if strcmp(algo.newton.solver, 'pcg')
    assert(algo.newton.system_opt == 4, "PCG can only be used for iterative system");
end

algo.newton.solver_used = algo.newton.solver;


algo.linesearch = struct;
algo.linesearch.time = 0;
algo.linesearch.total_time = 0;


%% initalize
algo.AL_penalty = params.AL_penalty_init;
if strcmp(params.method, 'drsip') || strcmp(params.method, 'fom')
    %% define functions
    % F = Z - X + model.c + Adagger * (A * (2 * X - Z - gamma * model.c) - model.b ) where gamma is the penalty parameter of augmented Lagrangian
    if model.n_box
        algo.funcF = @(x, z) z - x + algo.AL_penalty * model.c_int + algo.Adagger(algo.AXmap_int( x * 2 - z - algo.AL_penalty * model.c_int) -  model.b_int );
    else
        algo.funcF = @(x, z) z - x + algo.AL_penalty * model.c + algo.Adagger(algo.AXmap( x * 2 - z - algo.AL_penalty * model.c) -  model.b );
    end
    
    % proximal of h(x) = <c,x> + \delta_{Ax=b}
    % prox_h = (I - Adagger * A) (x - gamma * c) + Adagger * b
    algo.prox_h = @(x, gamma) (x - gamma * model.c) - algo.Adagger(algo.AXmap(x - gamma * model.c) - model.b);
end
if ~isempty(init_point)
    error("Not implemented yet");
else
    if strcmp(params.method, 'drsip')
        if params.init_strategy == 1
            [Z, mu, X, S, n_cones] = gen_init_point1_drspf(model, algo, params);
        elseif params.init_strategy == 2
            [Z, mu, X, S, n_cones] = gen_init_point2_drspf(model, algo, params);
        end
        if params.mu_init ~= -1
            mu = params.mu_init;
        elseif params.fom_start
            mu = 0;
        end
    elseif  strcmp(params.method, 'fom')
        if params.init_strategy == 1
            [Z, mu, X, S, n_cones] = gen_init_point1_drspf(model, algo, params);
        elseif params.init_strategy == 2
            [Z, mu, X, S, n_cones] = gen_init_point2_drspf(model, algo, params);
        end
        mu = 0;
    end
end

assert(VarSDP.isinterior(model.K, X));
assert(VarSDP.isinterior(model.K, S));

%% initialize history
hist = struct;
hist.pobjorg = [];
hist.dobjorg = [];
hist.pinforg = [];
hist.dinforg = [];
hist.gaporg = [];
hist.pvd = [];   % pinforg / dinforg
hist.dvp = [];   % dinforg / pinforg

%% start main loop
profile on
for iter = 0: params.max_iter
    algo.iter = iter;
    print_head = 0;
    
    %% decide method
    if strcmp(params.method, 'drsip')
        if params.fom_start > 0 && iter < params.fom_start
            algo.newton.method = 'fom' ;
            algo.newton.scaling = 'none';
        else
            algo.newton.method = 'drsip' ;
            algo.newton.scaling = 'none';
            if params.fom_start && iter == params.fom_start
                mu = 1e-3;
                fprintf(fid, "Switch to DRSPF\n");
                print_head = 1;
            end
        end
    elseif strcmp(params.method, 'fom')
        algo.newton.method = 'fom' ;
    else
        error("Invalid method: %s", params.method);
    end
    
    if iter == 0 || (params.fom_start && iter == params.fom_start)% initialize 
        [X, DX] = prox_neglog(Z, algo.AL_penalty * mu, model.K, params.det_reg, params.prox_neglog_formula);
        dZ = 0;
        stepsize = 1;
        n_failed_iter = 0;
    elseif params.AL_penalty_update_flag && mod(iter, params.AL_penalty_update_iter) == 0 % update AL penalty and do one first-order step
        AL_penalty_old = algo.AL_penalty;
        algo = AL_penalty_update(hist, algo, params);
        if algo.AL_penalty ~= AL_penalty_old
            %     rate_temp = algo.AL_penalty / AL_penalty_old;
            %     temp = rate_temp * (Z - X);
            %     Z = temp + algo.prox_h(X - temp, algo.AL_penalty);
            %     [X, V] = prox_neglog(Z, algo.AL_penalty * mu, model.K, params.det_reg);
            fprintf(fid, 'Augmented Lagrangian penalty is updated to %e\n', algo.AL_penalty);
            %     algo.funcF = @(x, z) z - x + algo.AL_penalty * model.c + algo.Adagger(algo.AXmap( x * 2 - z - algo.AL_penalty * model.c) -  model.b );
        end
        newt_flag = '1';
        error("Not implemented yet");
    elseif strcmp(algo.newton.method, 'fom') % do one first-order step
        Z = Z - F;
        [X] = prox_neglog(Z, algo.AL_penalty * mu, model.K, params.det_reg, params.prox_neglog_formula);
        

    elseif strcmp(algo.newton.method, 'drsip')
        %% compute r = - (F - mu * DmuF)
        beta = max([0.1, 1 - 0.1 * norm(F)]);
        beta = 1;
        S = X - Z; % S = algo.AL_penalty * mu * X^{-1}
        X0 = projection(Z, model.K);
        F0 = algo.funcF(X0, Z);
        temp = Jac_lmut(model.K, DX, S);
        muDmuF = - temp + algo.Adagger(algo.AXmap_int(temp * 2)); % muDmuF = mu * DmuF
        rhs = - ( F - beta * muDmuF );
        r = F - muDmuF - F0;
        Arhs = model.b_int - algo.AXmap_int(X - temp);
        % fprintf("norm(r) = %e, norm(F) = %e, norm(mu * DmuF) = %e\n", norm(r), norm(F), norm(mu * DmuF));
        
        if params.newton_reg
            algo.newton.sigma = algo.newton.kappa * (normF ^ algo.newton.sigPow);
        else
            algo.newton.sigma = 0;
        end
        
        [dZ, res_rd_pred, res_nt_pred, model, algo] = pred(X, DX, Z, mu, rhs, Arhs, model, algo, params);
        
        % res_augm = norm([DX(dZ) - dZ + algo.Atymap(w) + r ;algo.AXmap(DX(dZ) - r)] ) / (1 + norm([r; algo.AXmap(r)]));
        
        %% compute residue of the newton system, i.e.,
        % res = r - J * dZ where J is the Jacobian of F
        % temp = Jac_lmut(K, DX, dZ);
        % res_newt = norm(r - (dZ - temp + algo.Adagger(algo.AXmap(temp * 2 - dZ)))) / ( 1 + norm(r));
        

        %% linesearch and determine whether to accept the newton step
        % if res_nt_pred > 1e-2
        %     warning("newton system is solved badly");
        %     newt_flag = 'f';
        % else
            [stepsize, X_new, Z_new, model, algo] = linesearch_drsip_old(X, Z, dZ, mu, model, algo, params);
            if stepsize > 1e-5
                newt_flag = 's';
            else
                newt_flag =  'f';
            end
        % end

        %% update X, Z, mu
        if strcmp(newt_flag, 's')
            % newton step is successful, use second order method to update
            Z = Z_new;
            mu =  mu * (1 - beta * stepsize) ;
            
            %% another update
            n_failed_iter = 0;
        else
            % newton step is not successful, use first order method to update
            Z = Z - F;
            n_failed_iter = n_failed_iter + 1;
            algo.newton.method = 'fom' ;
        end
        
        [X, DX] = prox_neglog(Z, algo.AL_penalty * mu, model.K, params.det_reg, params.prox_neglog_formula);
        
        % %% update F
        % X0 = projection(Z, model.K);
        % F0 = algo.funcF(X0, Z);
        % X = prox_neglog(Z, algo.AL_penalty * mu, model.K, params.det_reg, params.prox_neglog_formula);
        % F = algo.funcF(X, Z);
        % r = F - muDmuF - F0;
        % %% update mu
        % if norm(F_new) <= max([algo.update.eta * phi, norm(F - F0) / algo.update.alpha1, norm(r) / algo.update.alpha2])
        %     mu = mu * algo.update.beta;

        % end
    end
    
    if params.accurate_recover
        X_prox = X;
        %% recover variables X, y, S using projection
        X = projection(Z, model.K);
        F = algo.funcF(X, Z);
        S = (X - Z) / algo.AL_penalty;
        y = invAAt(model.b - algo.AXmap(X + S - model.c), model, algo, params);
    else
        F = algo.funcF(X, Z);
        %% recover dual variables y, S
        S = (X - Z) / algo.AL_penalty;
        % y =  - (A * A')^{-1} * ( (A * X - b) / AL_penalty + A * (S - c) )
        if model.n_box
            y = invAAt(model.b_int / algo.AL_penalty - algo.AXmap_int(X / algo.AL_penalty + S - model.c_int), model, algo, params);
        else
            y = invAAt(model.b / algo.AL_penalty - algo.AXmap(X / algo.AL_penalty + S - model.c), model, algo, params);
        end
    end

        
    % update kappa
    normF_new = norm(F);
    if params.newton_reg && iter > 0
        % ratio =  - <F, dZ> / norm(dZ) ^ 2
        normdZ2 = norm(dZ)^2;
        if normF < 1e-3; normdZ2 = normF_new * normdZ2; end
        algo.newton.ratio = - inner_product(F, dZ) / normdZ2;
        
        if normF_new > normF
            algo.newton.kappa = algo.newton.kappa * 2;
            % continue;
        end
        
        %% update sigPow
        % if normF_new < 0.9 * normF
        algo.newton.sigPow = sigPow_update(normF);
    end
    normF = normF_new;



    % compute residual of box constraints
    boxres = 0;
    for p = 1 : length(model.K)
        cone = model.K{p};
        if strcmp(cone.type, 'b2l')
            boxres = boxres + sum(abs(X{p}(1: cone.size / 2) + X{p}(cone.size / 2 + 1 : end) - cone.params(:, 2) + cone.params(:, 1)));
        end
    end
    
    %% compute infeasibility of the current problen
    pres = algo.presmap_int(X);
    if model.n_box
        pres = pres(1: end - model.n_box);
    end
    pinfint = scale_b * norm(pres ) / (1 + scale_b * norm(model.b));
    dres =  algo.dresmap_int(y, S);
    % dres = 1 / algo.AL_penalty * (X - X_old);
    dinfint = scale_c * norm(dres) / ( 1 + scale_c * norm(model.c));
    if model.n_box
        pobjint = scale_b * scale_c * inner_product(X, model.c_int) + model.c0;
        y1 = y(1: end - model.n_box);
        dobjint = scale_b * scale_c * model.b' * y1 + model.c0; 
        for p = 1 : length(model.K)
            cone = model.K{p};
            if strcmp(cone.type, 'b2l')
                dobjint = dobjint -  scale_b * scale_c * S{p}(cone.size / 2 + 1 : end)' * (cone.params(:, 2) - cone.params(:, 1));
            end
        end
    else
        pobjint = scale_b * scale_c * inner_product(X, model.c) + model.c0;
        dobjint = scale_b * scale_c * model.b' * y + model.c0;
    end
    gapint = abs(pobjint - dobjint) / (1 + abs(pobjint) + abs(dobjint));
    % in theory , it holds that
    % pres = A * F
    % dres = - F / AL_penalty ;
    % pres + A * dres * AL_penalty = 0
    % gap = <A' * y + X, F> / AL_penalty + <X, S>
    % fprintf("pres - A * F = %e\n", norm(algo.AXmap_int(X) - model.b_int - algo.AXmap_int(F)));
    % fprintf("dres + F / AL_penalty = %e\n", norm(dres + F / algo.AL_penalty));
    % fprintf("pres / AL_penalty + A * dres = %e\n", norm(pres / algo.AL_penalty + algo.AXmap(dres)));
    % fprintf("gap - <A' * y + X / AL_penalty, F> - <X, S> = %e\n", abs(inner_product(X, model.c_int) - model.b_int' * y - inner_product(algo.Atymap_int(y) + X / algo.AL_penalty, F) - inner_product(X, S)));
    
    if params.accurate_recover
        X = X_prox; % remember to use the proximal solution for next iteration
    end

    % recover variables of original problem
    X_org = X;
    y_org = y;
    S_org = S;

    % recover from scale b, c
    X_org = scale_b * X_org;
    y_org = scale_c * y_org;
    S_org = scale_c * S_org;

    % recover from precond_A
    if params.precond_A
        y_org = algo.AAt.bwsolve(y_org) ;
    end

    %% recover variables genereated from b2l
    [X_org, y_org, S_org, dobjorg_box] = cone_b2l_recover(X_org, y_org, S_org, model_unscaled);

    %% recover variables genereated from u2l
    if params.u2l
        for p = 1: length(model.K)
            if strcmp(model.K{p}.type, 'u2l')
                cone = model.K{p};
                X_org{p} = X_org{p}(1: cone.size / 2) - X_org{p}(cone.size / 2 + 1: cone.size);
                S_org{p} = S_org{p}(1: cone.size / 2) ;
            end
        end
    end


    % recover from row_col_scale
    if params.rescale_option
        X_org = X_org .* algo.rescale.right_scale;
        y_org = y_org .* algo.rescale.left_scale;
        for p = 1: length(std_model.K)
            cone = std_model.K{p};
            if strcmp(cone.type, 'b')
                S_org{p}(1: cone.size) = S_org{p}(1: cone.size) ./ algo.rescale.right_scale{p};
                S_org{p}(cone.size + 1: end) = S_org{p}(cone.size + 1: end) ./ algo.rescale.right_scale{p};
            else
                S_org{p} = S_org{p} ./ algo.rescale.right_scale{p};
            end
        end
    end

    % recover from r2q
    [X_org, y_org, S_org] = cone_r2q_recover(X_org, y_org, S_org, model);
    
    
    % %% compute gaporg and infeasibility of the original problem
    pobjorg = inner_product(X_org, std_model.c);
    dobjorg = inner_product(std_model.b, y_org) + dobjorg_box;
    gaporg = abs(pobjorg - dobjorg) / (1 + abs(pobjorg) + abs(dobjorg));
    pinforg = norm(AXfun(std_model.K, std_model.At, X_org) - std_model.b) / (1 + norm(std_model.b));
    dinforg = norm(dresmap(std_model.K, std_model.At, y_org, S_org, std_model.c, [])) / (1 + norm(std_model.c));
    
    % pobjorg = pobjint;
    % dobjorg = dobjint;
    % gaporg = gapint;
    % pinforg = pinfint;
    % dinforg = dinfint;
    
    
    %% record history. since we start iter from 0, we need to add 1 to iter
    hist.pobjorg(iter+1) = pobjorg;
    hist.dobjorg(iter+1) = dobjorg;
    hist.gaporg(iter+1) = gaporg;
    hist.pinforg(iter+1) = pinforg;
    hist.dinforg(iter+1) = dinforg;
    hist.pvd(iter+1) = pinforg / dinforg;
    hist.dvp(iter+1) = dinforg / pinforg;
    
    if params.vector_mu
        mu_avg = MatCell.sum_all(mu) / n_cones;
    else
        mu_avg = mu;
    end
    
    %% print logs
    if params.print_log
        fm = '%.1e';
        fm_time = '%.1f';
        
        % basic log
        log_info = {'%6s', 'iter', '%6s', num2str(iter);
        % '%3s', 'flg', '%3s', newt_flag;
        '%6s', 'method', '%6s', algo.newton.method;
        '%7s', 'mu', '%7s', num2str(mu_avg, fm);
        '%8s', 'step', '%+8s', num2str(stepsize, fm);
        '%8s', 'pobj', '%+8s', num2str(pobjorg, fm);
        '%8s', 'dobj', '%+8s', num2str(dobjorg, fm);
        };

        % if true
        if params.precond_A || params.rescale_option
            log_info = [log_info;
                {
                '%8s', 'gap', '%+8s', num2str(gaporg, fm);
                '%8s', 'pinf', '%+8s', num2str(pinforg, fm);
                '%8s', 'dinf', '%+8s', num2str(dinforg, fm);
                % '%8s', 'pobjint', '%+8s', num2str(pobjint, fm);
                % '%8s', 'dobjint', '%+8s', num2str(dobjint, fm);
                % '%8s', 'gapint', '%+8s', num2str(gapint, fm);
                % '%8s', 'pinfint', '%+8s', num2str(pinfint, fm);
                % '%8s', 'dinfint', '%+8s', num2str(dinfint, fm);
                % '%8s', 'boxres', '%+8s', num2str(boxres, fm);
                }];
        else
            log_info = [log_info;
                {
                '%8s', 'gap', '%8s', num2str(gaporg, fm);
                '%8s', 'pinf', '%+8s', num2str(pinforg, fm);
                '%8s', 'dinf', '%+8s', num2str(dinforg, fm);
                }];
        end
        
        log_info = [log_info;
            {
            '%7s', 'time', '%7s', num2str(toc(t0), fm_time)
            }];
        
        % verbose log
        if params.print_log  == 2
            log_info = [log_info;
                {
                '%9s', 'solver', '%9s', algo.newton.solver_used;
                '%8s', 'dim', '%8d', algo.newton.dim;
                % '%8s', 'res_rd', '%+8s', num2str(res_rd_pred, fm);
                % '%8s', 'res_ag', '%+8s', num2str(res_augm, fm);
                % '%8s', 'res_pred', '%+8s', num2str(res_nt_pred, fm);
                '%8s', 'condH', '%+8s', num2str(algo.newton.condH, fm);
                }];
            if params.newton_reg > 0
                log_info = [log_info;
                    {
                    '%8s', 'sigma', '%+8s', num2str(algo.newton.sigma, fm);
                    '%8s', 'kappa', '%+8s', num2str(algo.newton.kappa, fm);
                    '%8s', 'sigPow', '%+8s', num2str(algo.newton.sigPow, fm);

                    }];
            end
            log_info = [log_info;
                {
                '%8s', 'X', '%+8s', num2str(norm(X), fm);
                '%8s', 'y', '%+8s', num2str(norm(y), fm);
                '%8s', 'S', '%+8s', num2str(norm(S), fm);
            %     % '%8s', 'XdotS', '%+8s', num2str(inner_product(X, S), fm);
                '%8s', 'Z', '%+8s', num2str(norm(Z), fm);
                '%8s', 'F', '%+8s', num2str(norm(F), fm);
                }];
            % if params.predcorr
            %     log_info = [log_info;
            %         {
            %         '%8s', 'dX_corr', '%+8s', num2str(norm(dX_corr), fm);
            %         '%8s', 'dS_corr', '%+8s', num2str(norm(dS_corr), fm);
            %         }];
            % end
            log_info = [log_info;
                {
                %                 '%8s', 'dim', '%8s', num2str(algo.newton.dim);
                '%7s', 'time_con', '%7s', num2str(algo.newton.time_construct, fm_time);
                '%7s', 'time_fac', '%7s', num2str(algo.newton.time_fact, fm_time);
                '%7s', 'time_slv', '%7s', num2str(algo.newton.time_solve, fm_time);
                % '%7s', 'time_ls', '%7s', num2str(algo.linesearch.time, fm_time)
                }];
            
            if strcmp(algo.newton.solver, 'pcg')
                log_info = [log_info;
                    {
                    '%5s', 'pcgit', '%5s', num2str(algo.newton.pcg_iter)
                    }];
            else
                log_info = [log_info;
                    {
                    % '%8s', 'rd_reg', '%8s', num2str(algo.newton.reg, fm) % regularization term added to the reduced system to make the cholesky or ldl factorizable
                    }];
            end
        end
        
        
        if mod(iter, 20) == 0 || print_head
            fprintf(fid, '\n');
            for i = 1:size(log_info, 1)
                fprintf(fid, log_info{i, 1}, log_info{i, 2});
                fprintf(fid, ' ');
            end
            fprintf(fid, '\n');
        end
        
        for i = 1:size(log_info, 1)
            fprintf(fid, log_info{i, 3}, log_info{i, 4});
            fprintf(fid, ' ');
        end
        fprintf(fid, '\n');
    end
    
    %% check status
    if strcmp(params.check_optimality_model, 'internal')
        acc= max([gapint, pinfint, dinfint]);
    elseif strcmp(params.check_optimality_model, 'original')
        acc = max([gaporg, pinforg, dinforg]);
    else
        error('unknown optimality model');
    end
    
    
    if acc < params.tol
        status = 'OPTIMAL';
    elseif n_failed_iter >= 20
        status = 'FAILED';
        %         elseif any(isnan(dZ)) || isnan(normF)
        %             status = 'NUMERICAL ERROR';
    elseif iter == params.max_iter
        status = 'MAX ITERATION REACHED';
    elseif toc(t0) > params.max_time
        status = 'TIME LIMIT EXCEEDED';
    else
        status = 'RUNNING';
    end
    
    
    if ~strcmp(status, 'RUNNING')
        break;
    end
end
profile off
%% print summary to log
out_summary = sprintf('%20s %20s\n%20s %20s\n%20s %20s\n%20s %20s\n%20s %20s\n%20s %20s\n%20s %20s\n', ...
    "Status", status, "Primal objective", num2str(pobjorg, '%.8e'), ...
    "Dual objective", num2str(dobjorg, '%.8e'), ...
    "Dual gap", num2str(gaporg, '%.2e'), ...de
    "Primal infeasibility", num2str(pinforg, '%.2e'), ...
    "Dual infeasibility", num2str(dinforg, '%.2e'), ...
    "Time (seconds)", num2str(toc(t0), '%.2f'),...
    "Iterations", num2str(iter));
fprintf(fid, '\n\n%s\n\n', out_summary);

%% close log file
if fid ~= 1
    fclose(fid);
end


if strcmp(algo.newton.solver, 'pardiso')
    pardisofree(algo.newton.pardiso_info);
end



%% print summary to screen
fid = 1;
fprintf(fid, '\n');
if isfield(model, 'name');   fprintf(fid, 'Problem            %s\n', model.name); end
fprintf(fid, 'Status             %s\n', status);
fprintf(fid, 'Primal Objective   %.8e\n', pobjorg);
fprintf(fid, 'Dual   Objective   %.8e\n', dobjorg);
fprintf(fid, 'Duality gap        %.2e\n', gaporg);
fprintf(fid, 'Primal infeas.     %.2e\n', pinforg);
fprintf(fid, 'Dual   infeas.     %.2e\n', dinforg);
fprintf(fid, 'Time               %.1f seconds\n', toc(t0));
fprintf(fid, 'Iterations         %d\n', iter);
if params.precond_A
    fprintf(fid, 'Precond A          %s\n', model.precond.solver);
end
fprintf(fid, '\n');

%% output the solution and other information

% out.Z = Z;
out.X = X_org;
out.y = y_org;
out.S = S_org;
out.mu = mu;

out.iter = iter;
out.pinf = pinforg;
out.dinf = dinforg;
out.gap = gaporg;
out.time = toc(t0);
out.status = status;

end






function params = set_default_params(params)
if ~isfield(params, 'max_iter'); params.max_iter = 100; end
if ~isfield(params, 'max_time'); params.max_time = 3600; end
if ~isfield(params, 'print_log'); params.print_log = 1; end % 0: no log, 1: standard log, 2: verbose log
if ~isfield(params, 'print_params'); params.print_params = 1; end % 1: print, 0: no print
if ~isfield(params, 'tol'); params.tol = 1e-6; end
if ~isfield(params, 'stepsize_init'); params.stepsize_init = 0.5; end
if ~isfield(params, 'warning_on'); params.warning_on = 1; end
if ~isfield(params, 'log_path'); params.log_path = ''; end
if ~isfield(params, 'newton_solver'); params.newton_solver = 'default'; end   % available solvers: 'ldl', 'ldlchol', 'chol', 'pcg', 'matlab_backslash', 'pardiso', 'default'
if ~isfield(params, 'AAt_solver'); params.AAt_solver = 'default'; end   % available solvers: 'ldlchol', 'chol', 'matlab_backslash', 'default'
% if ~isfield(params, 'mat_path'); params.mat_path = ''; end
% if ~isfield(params, 'solve_dual'); params.solve_dual = 0; end
% if ~isfield(params, 'solve_augmented'); params.solve_augmented = 0; end
if ~isfield(params, 'mu_init'); params.mu_init = -1; end % -1: auto 
% if ~isfield(params, 'path_following_coeff'); params.path_following_coeff = 1; end
if ~isfield(params, 'scale_bc_flag'); params.scale_bc_flag = 1; end
if ~isfield(params, 'scale_b'); params.scale_b = 0; end
if ~isfield(params, 'scale_c'); params.scale_c = 0; end
if ~isfield(params, 'precond_A'); params.precond_A = 0; end
if ~isfield(params, 'rescale_option'); params.rescale_option = 2; end
% if params.precond_A
%     params.rescale_option = 0;
%     fprintf('[warning] rescale_option is set to 0 because precond_A is set to 1\n');
% end
if ~isfield(params, 'system_opt'); params.system_opt = 0; end
if ~isfield(params, 'densemat_threshold'); params.densemat_threshold = 0.4; end
if ~isfield(params, 'newton_tol'); params.newton_tol = 1e-12; end
if ~isfield(params, 'newton_reg'); params.newton_reg = 0; end
if params.mu_init == 0
    if params.newton_reg == 0
        params.newton_reg = 1;
        fprintf('[warning] newton_reg is set to 1 because mu_init is set to 0\n');
    end
end
if ~isfield(params, 'AL_penalty_init'); params.AL_penalty_init = 1; end % penalty coeffient for augmented lagrangian function
if ~isfield(params, 'AL_penalty_update_flag'); params.AL_penalty_update_flag = 0; end
% if ~isfield(params, 'AL_penalty_update_iter'); params.AL_penalty_update_iter = 10; end
% if ~isfield(params, 'AL_penalty_min'); params.AL_penalty_min = 1e-6; end
% if ~isfield(params, 'AL_penalty_max'); params.AL_penalty_max = 1e6; end
% if ~isfield(params, 'AL_penalty_update_threshold'); params.AL_penalty_update_threshold = 1; end
% if ~isfield(params, 'AL_penalty_update_factor'); params.AL_penalty_update_factor = 1.3; end
if ~isfield(params, 'vector_mu'); params.vector_mu = 0; end
if ~isfield(params, 'check_optimality_model'); params.check_optimality_model = 'internal'; end
if ~isfield(params, 'scaling_method'); params.scaling_method = 'none'; end
if ~isfield(params, 'u2l'); params.u2l = 0; end
if ~isfield(params, 'method'); params.method = 'drsip'; end
assert(ismember(params.check_optimality_model, {'original', 'internal'}));
if ~isfield(params, 'det_reg'); params.det_reg = 0; end
if ~isfield(params, 'prox_neglog_formula'); params.prox_neglog_formula = 1; end
assert(ismember(params.prox_neglog_formula, [1, 2]));
if params.newton_reg
    if params.prox_neglog_formula == 2
        params.prox_neglog_formula = 1;
        warning('[warning] prox_neglog_formula is set to 1 because newton_reg is on');
    end
end
if ~isfield(params, 'init_strategy') 
    if strcmp(params.method, 'drsip') || strcmp(params.method, 'fom')
        params.init_strategy = 1; 
    else
        params.init_strategy = 3;
    end
end
if ~isfield(params, 'presolve'); params.presolve = 'none'; end
if ~isfield(params, 'fom_start'); params.fom_start = 0; end
if ~isfield(params, 'accurate_recover'); params.accurate_recover = 0; end
if ~isfield(params, 'detect_very_sparse_rows'); params.detect_very_sparse_rows = 0; end
if ~isfield(params, 'dense_column_strategy'); params.dense_column_strategy = 'default'; end

if params.system_opt == 2
    params.dense_column_strategy = 'none';
    params.detect_very_sparse_rows = 1;
    fprintf('[warning] dense_column_strategy is set to ''none'' and detect_very_sparse_rows is set to 1 because system_opt is set to 2\n');
end

end


function [AX, AXorg] = AXmap(X, K, At, precond)
m = size(At{1}, 2);
AX = zeros(m, 1);
for p = 1:length(K)
    cone = K{p};
    if isempty(At{p})
        continue;
    end
    Xp = X{p};
    if strcmp(cone.type, 's')
        % if isa(X{p}, "Var_sdp")  % X{p} is Var_sdp object
        %    x = mysvec(X{p});
        %    AX = AX + (x'*At{p})';
        % else % X{p} is one single matrix or vector.
        if length(cone.size) > 1
            Xp = sparse(Xp);
        end
        if size(X{p}, 1) == sum(cone.size .* (cone.size + 1) / 2) % X{p} is in vector form
            AX = AX + At{p}' * Xp;
        elseif size(X{p}, 1) == sum(cone.size) % X{p} is in matrix form
            x = mysvec(cone, Xp);
            AX = AX + At{p}' * x;
        else
            error('AXfun: dimension mismatch for block %i', p);
        end
    elseif strcmp(cone.type, 'b2l')
        assert(mod(length(Xp), 2) == 0);
        AX = AX + At{p}' * Xp(1: length(Xp) / 2);
    else
        AX = AX + At{p}' * Xp;
    end
end
AXorg = AX;
if isfield(precond, 'fwsolve')
    AX = precond.fwsolve(AXorg);
else
    AX = fwsolve(precond, AXorg);
end
end

function Aty = Atymap(y, K, At, precond)
if isfield(precond, 'bwsolve')
    Aty = Atyfun(K, At, precond.bwsolve(y));
else
    Aty = Atyfun(K, At, bwsolve(precond, y));
end
end

function [pres] = presmap(K, At, X, b, precond)
pres = AXmap(X, K, At, precond) - b;
end

function [dres] = dresmap(K, At, y, S, c, precond)
if isempty(precond)
    dres = Atyfun(K, At, y);
else
    dres = Atymap(y, K, At, precond);
end

for p = 1: length(K)
    cone = K{p};
    if strcmp(cone.type, 'b2l') 
        dres{p} = dres{p} + S{p}(1: cone.size / 2) - S{p}(cone.size / 2 + 1: cone.size) - c{p};
    elseif strcmp(cone.type, 'b') 
        dres{p} = dres{p} + S{p}(1: cone.size) - S{p}(cone.size + 1: cone.size * 2)- c{p};
    else
        dres{p} = dres{p} + S{p} - c{p};
    end
end
end

function [AX_int, AX, res_box] = AXmap_int(X, K, At, precond)
AX = AXmap(X, K, At, precond);
AX_int = AX;
res_box = [];
for p = 1: length(K)
    cone = K{p};
    if strcmp(cone.type, 'b2l')
        res_box = [res_box;
            X{p}(1: cone.size / 2) + X{p}(cone.size / 2 + 1: cone.size)];
    end
end
AX_int = [AX_int; res_box];
end

function [pres] = presmap_int(K, At, X, b, precond)
pres = AXmap(X, K, At, precond) - b;
for p = 1: length(K)
    cone = K{p};
    if strcmp(cone.type, 'b2l')
        pres = [pres;
            X{p}(1: cone.size / 2) + X{p}(cone.size / 2 + 1: cone.size) - cone.params(:, 2) + cone.params(:, 1)];
    end
end
end

function [Aty] = Atymap_int(K, At, y, precond, n_box)
    if nargin < 5
        n_box = 0;
    end
    y_box = y(end - n_box + 1: end);
    y = y(1: end - n_box);
    Aty = Atymap(y, K, At, precond);
    if n_box == 0
        return;
    end
    idx = 1;
    for p = 1: length(K)
        cone = K{p};
        len = cone.size / 2;
        if strcmp(cone.type, 'b2l')
            Aty{p} = [Aty{p} + y_box(idx: idx + len - 1);
                      y_box(idx: idx + len - 1)];
        else
            
        end
        idx = idx + len;
    end
end

