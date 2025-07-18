%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-06 13:19:08
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
%% socp + lp
%% min <c_q, x_q> + <c_l, x_l>
% st    x_q : second order cone
%       x_l >= 0
%       A_q * x_q + A_l * x_l == b

clear
randn('seed', 1);
rand('seed',1); 
m = 1000;
n_q = 2000;
n = n_q ;


A_q = randn(m, n_q);
% A_q = sprand(m, n_q, 0.4);
c_q = rand(n_q, 1);
c_q(1) = norm(c_q) + abs(rand(1)); % ensure the problem is bounded
x_q0 = rand(n_q ,1);
b_q = A_q * x_q0;



c = [c_q];
A = [A_q];
% assert(rank(A) == min([m, n]));
x0 = [x_q0];
b = b_q ;


% norm(A_sc * x - b_sc) <= d_sc' * x + gamma_sc
% 
% A_sc = diag([0; ones(n_q-1, 1)]);
% b_sc = zeros(n ,1);
% d_sc = [1; zeros(n-1, 1)];
% gamma_sc = 0;
% socConstraints = secondordercone(A_sc,b_sc,d_sc,gamma_sc);
% 
% fprintf("Calling matlab coneprog\n");
% [x, fval] = coneprog(c,socConstraints, ...
%                     [],[], ... % inequality condtraints A_ineq * x  <= b_ineq
%                     A,b,   ... % equation constrains A_eq * x == b_eq
%                     [0; -inf*ones(n_q-1, 1)], inf*ones(n, 1)) % lower bound and upper bound of x
% assert(x(1) >= norm(x(2:n_q)));

fprintf(repmat('*', 40, 1));
fprintf('\n');

% At = [A_q'; A_l'];
% At = sparse(At);
% C = [c_q; c_l];
% blk = {'q', n_q;
%        'l', n_l};

model.At = MatCell({A_q'});
model.c = MatCell({c_q});
model.b = b;
model.K = Cone({BasicCone('q', n_q)});
model.name = "random QP" ;
params.warning_on = 0;
% params.print_log = 2;
% params.newton_reg = 1;
params.max_iter = 200;
% params.system_opt = 3;
% params.predcorr = 0;
% params.scaling_method = "HKM" ;
% params.newton_reg = 1;
% params.rescale_option = 2;
params.prox_neglog_formula = 1;
% params.dense_column_strategy = 'augmented';
params.precond_A = 1;
% params.newton_solver='pardiso';
% params.AL_penalty_init = 1;
out  = drspf(model, params);

% addpath_sdpt3
% options = struct;
% options.vers = 1; % 1: HKM, 2: NT
% options.predcorr = 0;
% options.printlevel=inf;
% out2 = sqlp(Cone.toblk(model.K), model.At, model.c, model.b, options);


