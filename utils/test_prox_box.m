%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-12 21:00:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-12 21:45:02
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
mu = 10;

n = 1e5;
l = randn(n, 1);
u = l + 10 * abs(rand(n, 1));
z = randn(n, 1);
tic;
x = prox_box(z, mu, l, u, 3);
toc;

%% check by x - z - \mu / ( x - l) + \mu / (u - x) = 0
res = x - z - mu ./ (x - l) + mu ./ (u - x);
fprintf('check by x - z - mu / ( x - l) + mu / (u - x) = 0, error = %e\n', norm(res));


% tic;
% x = prox_box(z, mu, l, u, 2);
% toc;
% 
% %% check by x - z - \mu / ( x - l) + \mu / (u - x) = 0
% res = x - z - mu ./ (x - l) + mu ./ (u - x);
% fprintf('check by x - z - mu / ( x - l) + mu / (u - x) = 0, error = %e\n', norm(res));


