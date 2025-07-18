%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-01-03 10:49:57
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
%% socp + lp
%% min <c_q, x_q> + <c_l, x_l>
% st    x_q : second order cone
%       x_l >= 0
%       A_q * x_q + A_l * x_l == b

rand('seed',1); 
m = 30;
n_q = 50;
n_l = 40;
n = n_q + n_l;

x_q0 = rand(n_q ,1);
x_q0(1) = norm(x_q0); % ensure abs(x0(1)) >= norm(x(2:end));
A_q = rand(m, n_q);
c_q = rand(n_q, 1);
c_q(1) = norm(c_q); % ensure the problem is bounded
b_q = A_q * x_q0;


x_l0 = rand(n_l ,1);
A_l = rand(m, n_l);
c_l = rand(n_l, 1);
b_l = A_l * x_l0;

c = [c_q; c_l];
A = [A_q, A_l];
x0 = [x_q0; x_l0];
b = b_q + b_l;


% % norm(A_sc * x - b_sc) <= d_sc' * x + gamma_sc
% A_sc = diag([0; ones(n_q-1, 1); zeros(n_l, 1)]);
% b_sc = zeros(n ,1);
% d_sc = [1; zeros(n-1, 1)];
% gamma_sc = 0;
% socConstraints = secondordercone(A_sc,b_sc,d_sc,gamma_sc);

% [x, fval] = coneprog(c,socConstraints, ...
%                     [],[], ... % inequality condtraints A_ineq * x  <= b_ineq
%                     A,b,   ... % equation constrains A_eq * x == b_eq
%                     [0; -inf*ones(n_q-1, 1); zeros(n_l, 1)], inf*ones(n, 1)) % lower bound and upper bound of x
% assert(x(1) >= norm(x(2:n_q)));


fprintf(repmat('*', 40, 1));
fprintf('\n');

At = [A_q'; A_l'];
At = sparse(At);
C = [c_q; c_l];
blk = {'q', n_q;
       'l', n_l};

model.At = MatCell({A_q', A_l'});
model.c = MatCell({c_q, c_l});
model.b = b;
model.K = Cone({BasicCone('q', n_q), BasicCone('l', n_l)});
model.name = "random QLP" ;
params.warning_on = 0;
params.print_log = 2;
params.init_stepsize = 1;
params.newton_reg = 0;
out  = drsip(model, params);


