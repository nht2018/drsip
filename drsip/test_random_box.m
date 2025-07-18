%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-10 16:38:55
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-19 11:56:37
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
%% lp with box
%% min <c, x>
% st    lb <= x <= ub
%       A * x == b

rand('seed',1); 
clear; 
m = 30;
n_l = 50;
n_b = 50;
n = n_l + n_b;



x = rand(n_l + n_b, 1);
A = sprand(m, n_l + n_b, 0.5);
% 
% % generate A such that A * A' = I
% B = rand(m, n);
% [~, ~, V] = svd(B);
% A = V(:, 1:m)';

c = rand(n, 1);
b = A * x;
x1 = x(1: n_b) ;
lb = x1 - 2 * abs(rand(n_b, 1)) - 1;
ub = x1 + 2 * abs(rand(n_b, 1)) + 1;



if n_l > 0 && n_b> 0
    model = struct();
    model.At = MatCell({ A(:, 1:n_b)', A(:, n_b+1:end)' });
    model.c = MatCell({c(1:n_b), c(n_b+1:end)});
    model.b = b;
    model.K = Cone({BasicCone('b', n_b, [lb, ub]), BasicCone('l', n_l)});
elseif n_l == 0
    model = struct();
    model.At = MatCell({ A' });
    model.c = MatCell({c});
    model.b = b;
    model.K = Cone({BasicCone('b', n_b, [lb, ub])});
else
    model = struct();
    model.At = MatCell({ A' });
    model.c = MatCell({c});
    model.b = b;
    model.K = Cone({BasicCone('l', n_l)});
end


model.name = "random box LP";
params = struct();
% params.newton_solver = 'chol';
% params.system_opt = 3;
params.max_iter = 200;
params.rescale_option = 0;
params.scale_bc_flag = 1;
params.scale_b = 1;
params.scale_c = 10;
params.precond_A = 0;
% params.AAt_solver = 'ldlchol';
% params.method = 'fom';
out = drspf(model, params);


% call gurobi
lb = [lb; zeros(n_l, 1)];
ub = [ub; inf(n_l, 1)];
[At, b, c, lb, ub] = row_col_rescale_mat(A', b, c, lb, ub,2);
A = At';

addpath_gurobi;
model_gurobi = struct();
model_gurobi.A = A;
model_gurobi.obj = c;
model_gurobi.rhs = b;
model_gurobi.lb = lb;
model_gurobi.ub = ub;
model_gurobi.sense = '=';
model_gurobi.vtype = 'C';
model_gurobi.modelsense = 'min';

out = gurobi(model_gurobi);
