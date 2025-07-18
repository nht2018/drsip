%%*****************************************************************
%% drspred: solve the Newton equation and obtain the direction
%%*****************************************************************

function [dZ, res_rd, res_nt, model, algo] = pred(X, DX, Z, mu, r, Ar, model, algo, params)

assert(length(r) == length(Z), "r length is not correct");

At = model.At;
K = model.K;


%% system_opt
% 1: ldl on sparse matrix
% 2: ldl on sparse matrix with mat11 PSD(deprecated)
% 3: chol on dense matrix
% 4: iterative method
% 5: augmented system(deprecated)



sigma = algo.newton.sigma; % regularization parameter for newton system

m = size(At{1}, 2);

t = tic;
%% initialize lhs
lhs = struct();


if ismember(algo.newton.system_opt, [1, 3]) % direct method
    lhs.mat11 = sparse(m, m);
    if ismember(algo.newton.system_opt, [1, 3])
        lhs.mat21 = sparse(0, m);
    end
    lhs.mat31 = sparse(0, m);
    lhs.mat32 = []; % lhs.mat32 is always empty
    lhs.mat22_diag = []; % handle free variables
    lhs.mat33 = []; % handle dense column
elseif algo.newton.system_opt == 4 % iterative method
    lhs.mat11 = MatCell(length(Z)); % for iterative method, mat11 is a cell storing each block
end

if ~isfield(algo.model, 'idx_u') 
    algo.model.idx_u = []; % block index of unbounded cone
    algo.model.size_u = [];
    for p = 1: length(model.K)
        cone = model.K{p};
        if strcmp(cone.type, 'u')
            algo.model.idx_u = [algo.model.idx_u, p];
            algo.model.size_u = [algo.model.size_u, cone.size];
        end
    end
end

%% compute H and construct lhs
algo.newton.q_dense_col = 0;

%% we find that divide 0 may occur in the computation of H, so we need to handle it in the future
% Htemp = (DX - (1+sigma)*I)^{-1}
% H = - (1+2*sigma) * (DX - (1+sigma)*I)^{-1} - I
if sigma > 0 || params.prox_neglog_formula == 1
    % if true
    Htemp = Jac_ops(DX, K, 'affine_inv', 1, -(1+sigma)); % this is also used for recover dZ from w
    H = Jac_ops(Htemp, K, 'affine', -(1+2*sigma), -1);
    for p = 1: length(K)
        if strcmp(K{p}.type, 'u')
            H{p}.shift = 0; % actually H{p}.shift= inf, since we do no need to compute H for unbounded cone, we set it to 0 for convenience
        end
    end
else % for sigma = 0, we cans compute Htemp = (DX - I)^{-1} and H = - Htemp - I = directly
    H = StructCell(length(K));
    for p = 1:length(K)
        cone = K{p};
        if params.vector_mu
            mup = mu{p} * algo.AL_penalty;
        else
            mup = mu * algo.AL_penalty;
        end
        if strcmp(cone.type,'l') || endsWith(cone.type,'l')
            H{p}.shift = X{p} .^2 ./ mup;
        elseif strcmp(cone.type, 's')
            assert(algo.newton.system_opt == 4)
            H{p}.lmut = @(r_) 1 / mup * X{p} * r_ * X{p}; % here r_ is a matrix
        elseif strcmp(cone.type,'q') || endsWith(cone.type,'q')
            detx = soc_ops(X{p}, 'det', cone.size);
            detx = max(detx, 1e-16); % avoid divide 0
            ind_head = cumsum(cone.size) - cone.size + 1;  % record the start index of each cone blocks
            tempdiag = ones(size(Z{p}));
            tempdiag(ind_head) = -1;
            H{p}.shift = (repelem( detx ./ mup, cone.size, 1)) .* tempdiag ;
            H{p}.coeff = detx ./ mup;
            H{p}.lr = sqrt( 2 ./ repelem(detx, cone.size, 1)) .* X{p};
        elseif strcmp(cone.type,'u')
            H{p}.shift = 0;
        end
    end

end
% H^{-1} = (1+2*sigma) * (DX + sigma * I)^{-1} - I
maxH = 0;
minH = inf;
for p = 1:length(K)
    if ~ strcmp(K{p}.type, 'u') && ~ strcmp(K{p}.type, 's') % case of s to be implemented
        maxH = max(maxH, max(abs(H{p}.shift)));
        minH = min(minH, min(abs(H{p}.shift)));
    end
end

algo.newton.condH = maxH / minH;



%% handle box constraints
exist_box = 0;
H_rd = H;
H1 = [];
H2 = [];
for p = 1:length(K)
    cone = K{p};
    if strcmp(cone.type, 'b2l')
        exist_box = 1;
        % reduce the dim from cone.size to cone.size / 2
        h1 = H{p}.shift(1: cone.size / 2);
        h2 = H{p}.shift(cone.size / 2 + 1: cone.size);
        H_rd{p}.shift = h1 .* h2 ./ (h1 + h2);
        H1 = [H1; h1];
        H2 = [H2; h2];
    end
end
n_box = model.n_box;
At_box = model.At_box;

[lhs, model, algo] = constructlhs(H_rd, model, algo, params);
% if algo.iter == 1
%     fprintf(algo.fid, "[newton system] dim1 = %d, dim2 = %d, dim3 = %d\n", lhs.dim1, lhs.dim2, lhs.dim3);
%     if algo.newton.system_opt == 1 || algo.newton.system_opt == 3
%         fprintf(algo.fid, "\tnonzeros = %d, sparsity = %f\n", nnz(lhs.mat), nnz(lhs.mat) / numel(lhs.mat));
%     end
% end
algo.newton.time_construct = toc(t);


%% factorize the reduced system
t = tic;
algo.newton.ordering = [];
[algo.newton] = linsysfact(lhs, algo.newton);
algo.newton.time_fact = toc(t);
algo.newton.total_time_fact = algo.newton.total_time_fact + algo.newton.time_fact;


t = tic;
%% compute r in the reduced system
% [~, rhs_A, rhs_box] = algo.AXmap_int(r - Jac_lmut(K, H, r));
r_temp = AXfun(K, model.At_int, r - Jac_lmut(K, H, r));
% r_temp = Ar - AXfun(K, model.At_int, Jac_lmut(K, H, r));
rhs_A = r_temp(1: end - n_box);
rhs_box = r_temp(end - n_box + 1: end);
if exist_box
    rhs_A = rhs_A - At_box' * spdiag( H1 ./ (H1 + H2)) * rhs_box;
end

if algo.newton.system_opt == 1
    rhs_u = MatCell.vert_concat(r, algo.model.idx_u);    % concatenation of dense part and unbouded blocks
    if strcmp(algo.newton.dense_column_strategy, 'augmented')
        rhs = [rhs_A;
            - rhs_u;
            zeros(lhs.dim3, 1)];
    else % SMW
        rhs = [rhs_A;
        - rhs_u];
    end
else % algo.newton.system_opt == 3 || algo.newton.system_opt == 4
    rhs_u = MatCell.vert_concat(r, algo.model.idx_u);  % rhs of unbounded blocks
    rhs = [rhs_A;
        - rhs_u];
end

%% solve the reduced system
[d, res_rd, algo] = solve_reduced_system(lhs, rhs, model, algo, params) ;

%% recover variables in Newton system
%% compute w
if algo.newton.system_opt ~= 2
    w = d(1: lhs.dim1); % if preconditioner is used, then d(1: lhs.dim1) actually is R^{-1} * w
else
    w = d;
end

d_u = d(lhs.dim1 + 1: lhs.dim1 + lhs.dim2);

if exist_box
    d_box = 1 ./ (H1 + H2) .* (rhs_box - H1 .* (At_box * w));
    w = [w; d_box];
else
    d_box = [];
end


%% compute dZ
temp = Atyfun(K, model.At_int, w) + r;
if sigma > 0
    dZ = - Jac_lmut(K, Htemp, temp);
else
    dZ = temp + Jac_lmut(K, H, temp);
end

if numel(algo.model.idx_u) > 1
    dZ{algo.model.idx_u} = MatCell.vert_split(d_u / (1+sigma), algo.model.size_u) ;
elseif numel(algo.model.idx_u) == 1
    dZ{algo.model.idx_u} = d_u / (1 + sigma);
end




%% check residual of Newton system
temp = Jac_lmut(K, DX, dZ);
res_nt = norm(r - (dZ - temp + algo.Adagger(algo.AXmap_int(temp * 2 - dZ)))) / ( 1 + norm(r));
% res_augm = norm([temp - dZ + algo.Atymap_int(w) + r ;
%                  algo.AXmap_int(temp - r)           ] ) / (1 + norm([r; algo.AXmap_int(r)])) ;




algo.newton.time_solve = toc(t);
algo.newton.total_time_solve = algo.newton.total_time_solve + algo.newton.time_solve;


end