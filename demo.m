startup;
rand('seed',1); 
m = 30; n_q = 50; n_l = 40; n = n_q + n_l;

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

model.At = MatCell({A_q', A_l'});
model.c = MatCell({c_q, c_l});
model.b = b_q + b_l;
model.K = Cone({BasicCone('q', n_q), BasicCone('l', n_l)});
model.name = "random QLP" ;
params = struct;
out  = drsip(model, params);