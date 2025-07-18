% Generate portfolio optimization problem
% max \mu^T * z -\gamma z^T * \Sigma *z 
% s.t. \mathbf{1}^T z=1,  z >= 0
% with \Sigma = F * F^T + D

run_drsip = 1;
run_ecos = 1;
run_scs = 1;
run_superscs = 1;
run_sdpnalplus = 1;
run_sdpt3 = 1;

% n = 10000;
% m = 2000;

if strcmp(computer, "MACI64") || strcmp(computer, "MACA64")
    log_path = sprintf("../logs/portfolio/n%d_m%d", n, m);
    result_path = sprintf("../results/portfolio/n%d_m%d", n, m);
elseif strcmp(computer, 'GLNXA64') % my linux server
    log_path = sprintf("../logs_linux/portfolio/n%d_m%d", n, m);
    result_path = sprintf("../results_linux/portfolio/n%d_m%d", n, m);
end

if ~exist(log_path, 'dir')
    mkdir(log_path);
end

if ~exist(result_path, 'dir')
    mkdir(result_path);
end


rng('default');
rng(1);
tol = 1e-6;

density = 0.1;
rc = 0.5; % estimated reciprocal condition number



mu = exp(0.01 * randn(n, 1)) - 1; % returns
D = rand(n,1) / 10; % idiosyncratic risk
F = sprandn(n, m, density, rc) / 10; % factor model
gamma = 1;


% general model with x = [z; u; v; t; s]
gen_model = struct;
gen_model.c0 = 0;
gen_model.c = [-mu; 0; 0; gamma; gamma];
gen_model.A = [ones(1, n), 0, 0, 0, 0];
gen_model.lc = 1;
gen_model.uc = 1;
gen_model.lx = [zeros(n, 1); -inf; -inf; -inf; -inf];
gen_model.ux = [inf * ones(n, 1); inf; inf; inf; inf];
gen_model.F = [sparse(1, n), 1, 0, 0, 0;
                spdiag(D .^ 0.5), sparse(n, 4);
                sparse(1, n), 0, 1, 0, 0;
                F', sparse(m, 4);
                sparse(1, n), 0, 0, 1, 0;
                sparse(1, n), 0, 0, -1, 0;
                sparse(1, n), 2, 0, 0, 0;
                sparse(1, n), 0, 0, 0, 1;
                sparse(1, n), 0, 0, 0, -1;
                sparse(1, n), 0, 2, 0, 0];
gen_model.g = [zeros(n+1, 1);
                zeros(m+1, 1);
                [1; 1; 0];
                [1; 1; 0]];
gen_model.D = {BasicCone('q', [n+1, m+1, 3, 3])};
gen_model.name = "portfolio" ;

std_model = standardize_forward(gen_model);



if run_drsip
params.warning_on = 0;
params.print_log = 2;
params.stepsize_init = 0.85;
params.rescale_option = 2;
params.newton_reg = 0;
params.tol = tol;
% out  = drsip(gen_model, params);
try
diary(fullfile(log_path, "drsip.log"));
tic;
out = drsip(std_model, params);
t = toc;
diary off;
if strcmp(out.status, 'OPTIMAL')
    status = 's';
else
    status = 'f';
end
fileID = fopen(fullfile(result_path, "drsip.txt"), 'w');
fprintf(fileID, '%.4e\t%i\t%.4e\t%.4e\t%.4e\t%s', t, out.iter, out.gap, out.pinf, out.dinf, status);
fclose(fileID);
catch
fileID = fopen(fullfile(result_path, "drsip.txt"), 'w');
fprintf(fileID, '%.4e\t%i\t%.4e\t%.4e\t%.4e\t%s', -1, -1, -1, -1, -1, 'f');
fclose(fileID);
end
end


if run_ecos
addpath_ecos;
% [x,y,info,s,z] = ecos(c,G,h,dims,A,b) Solves a pair of primal and 
% dual cone programs. The primal problem is defined as
%      minimize    c'*x
%      subject to  G*x <=_K h
%                  A*x = b
%    where c,G,h,dims are defined as above, and A is a sparse matrix of 
%    size (p,n), and b is a dense matrix of size (p,1).
%    It is assumed that rank(A) = p and rank([A; G]) = n.
G = [- speye(n), sparse(n, 4);
    - gen_model.F];
h = [zeros(n, 1);
    gen_model.g];
dims = struct;
dims.l = n;
dims.q = [n+1, m+1, 3, 3];
c = gen_model.c;
A = sparse(gen_model.A);
b = gen_model.lc;

% Optimize the problem. 
ecosopts = ecosoptimset('reltol',tol,'maxit',500);
try
diary(fullfile(log_path, "ecos.log"));
tic;
[x, y, info,s,z] = ecos(c, G, h, dims, A, b, ecosopts);
t = toc;
diary off;
if strcmp(info.infostring, 'Optimal solution found')
    status = 's';
else
    status = 'f';
end
fileID = fopen(fullfile(result_path, "ecos.txt"), 'w');
fprintf(fileID, '%.4e\t%i\t%.4e\t%.4e\t%.4e\t%s', t, info.iter, info.relgap, info.pres, info.dres, status);
fclose(fileID);
catch
fileID = fopen(fullfile(result_path, "ecos.txt"), 'w');
fprintf(fileID, '%.4e\t%i\t%.4e\t%.4e\t%.4e\t%s', -1, -1, -1, -1, -1, 'f');
fclose(fileID);
end


% addpath_cvx;
% cvx_begin
% cvx_solver ecos
% cvx_solver_settings('reltol', tol, 'feastol', tol, 'abstol', tol)
% variable x(n)
% maximize (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))
% sum(x) == 1
% x >= 0
% cvx_end


% x_viol = min(x)
% budget_viol = abs(1-sum(x))
% obj = (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)));


end


if run_scs
addpath_cvx;
addpath_scs;
cvx_setup
try
diary(fullfile(log_path, "scs.log"));
tic;
cvx_begin
cvx_solver scs
cvx_solver_settings('eps_rel',1e-6,'eps_infeas', 1e-6, 'scale',1)
variable x(n)
maximize (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))
sum(x) == 1
x >= 0
cvx_end
t = toc;
diary off;

x_viol = min(x);
budget_viol = abs(1-sum(x));
obj = (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)));
if cvx_status == "Solved"
    status = 's';
else
    status = 'f';
end
fileID = fopen(fullfile(result_path, "scs.txt"), 'w');
fprintf(fileID, '%.4e\t%i\t%.4e\t%.4e\t%.4e\t%s', t, cvx_slvitr, -1, -1, -1, status);
fclose(fileID);
catch
fileID = fopen(fullfile(result_path, "scs.txt"), 'w');
fprintf(fileID, '%.4e\t%i\t%.4e\t%.4e\t%.4e\t%s', -1, -1, -1, -1, -1, 'f');
fclose(fileID);
end
end

if run_superscs
addpath_cvx;
addpath_superscs;
cvx_setup
try
diary(fullfile(log_path, "superscs.log"));
tic;
cvx_begin
cvx_solver scs
cvx_solver_settings('eps', 1e-8, 'do_super_scs', 1, 'rho_x', 1,...
        'direction', 100, 'memory', 50);
variable x(n)
maximize (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))
sum(x) == 1
x >= 0
cvx_end
t = toc;
diary off;

x_viol = min(x);
budget_viol = abs(1-sum(x));
obj = (mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)));
if cvx_status == "Solved"
    status = 's';
else
    status = 'f';
end
fileID = fopen(fullfile(result_path, "superscs.txt"), 'w');
fprintf(fileID, '%.4e\t%i\t%.4e\t%.4e\t%.4e\t%s', t, cvx_slvitr, -1, -1, -1, status);
fclose(fileID);
catch
fileID = fopen(fullfile(result_path, "superscs.txt"), 'w');
fprintf(fileID, '%.4e\t%i\t%.4e\t%.4e\t%.4e\t%s', -1, -1, -1, -1, -1, 'f');
fclose(fileID);
end
end


if run_sdpnalplus
addpath_sdpnalplus;
blk = Cone.toblk(std_model.K);
try
diary(fullfile(log_path, "sdpnalplus.log"));
options = struct;
options.gaptol = tol;
tic;
[obj,X,y,S,~,~,~,~,info,runhist] = sdpnalplus(blk, std_model.At, std_model.c, std_model.b, -inf, inf, [], [], [], options);
t = toc;
diary off;
removepath_sdpnalplus;
if info.termcode == 0
    status = 's';
else
    status = 'f';
end
fileID = fopen(fullfile(result_path, "sdpnalplus.txt"), 'w');
fprintf(fileID, '%.4e\t%i\t%.4e\t%.4e\t%.4e\t%s', t, info.iter, info.relgap, runhist.primfeasorg(1), runhist.dualfeasorg(1), status);
fclose(fileID);
catch
fileID = fopen(fullfile(result_path, "sdpnalplus.txt"), 'w');
fprintf(fileID, '%.4e\t%i\t%.4e\t%.4e\t%.4e\t%s', -1, -1, -1, -1, -1, 'f');
fclose(fileID);
end
end

if run_sdpt3
addpath_sdpt3;
blk = Cone.toblk(std_model.K);
try
diary(fullfile(log_path, "sdpt3.log"));
options = struct;
options.gaptol = tol;
[obj,X,y,S,info,runhist] = sdpt3(blk, std_model.At, std_model.c, std_model.b, options);
diary off;
removepath_sdpt3;
if info.termcode == 0
    status = 's';
else
    status = 'f';
end
fileID = fopen(fullfile(result_path, "sdpt3.txt"), 'w');
fprintf(fileID, '%.4e\t%i\t%.4e\t%.4e\t%.4e\t%s', t, info.iter, info.relgap, info.pinfeas, info.dinfeas, status);
fclose(fileID);
catch
fileID = fopen(fullfile(result_path, "sdpt3.txt"), 'w');
fprintf(fileID, '%.4e\t%i\t%.4e\t%.4e\t%.4e\t%s', -1, -1, -1, -1, -1, 'f');
fclose(fileID);
end
end