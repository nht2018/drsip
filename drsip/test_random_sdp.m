%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-01-26 21:55:44
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-23 09:36:42
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
%% generate a random SDP problem
rng('default')
% rng(1); % warning: not all rand seed lead to sucuessful solution

K = Cone( { 
          BasicCone('s', [3 3 3]), ...
%           BasicCone('q', [5]), ...
          BasicCone('l', [5]) 
           });

m = 3;


At = MatCell(length(K));
c = MatCell(length(K));
x0 = MatCell(length(K));
for p=1: length(K)
    cone = K{p};
    if strcmp(cone.type, 's')
        n = sum(cone.size .* (cone.size + 1) / 2);
        At{p} = zeros(n, m);
        for i=1: m+2
            x = cell(length(cone.size), 1);
            for j=1: length(cone.size)
                x{j} = sparse(rand_psd(cone.size(j)));
            end
            if i <= m
                At{p}(:, i) = mexsvec({'s', cone.size}, blkdiag(x{:}), 1);
            elseif i == m+1
                c{p} = blkdiag(x{:});
            else
                x0{p} = blkdiag(x{:});
            end
        end

    else
        n = sum(cone.size);
        At{p} = randn(n, m);
        c{p} = randn(n, 1);
        x0{p} = randn(n, 1);
        if strcmp(cone.type, 'q')
            c{p}(1) = norm(c{p});  
            x0{p}(1) = norm(x0{p}); % make sure the problem is feasible        
        end
    end

end


b = AXfun(K, At, x0);

model = struct();
model.K = K;
model.At = At;
model.c = c;
model.b = b;

tol = 1e-6;
% 
% %% solve the problem by sdpt3
% addpath_sdpt3;
% options = struct;
% options.gaptol = tol;
% options.inftol = tol;
% [obj,X,y,Z,algo,runhist] = sqlp(Cone.toblk(K), At.data, c.data, b, options);


%% solve the problem
params.rescale_option = 0;
params.tol = tol;
params.print_log = 2;
params.max_iter = 50;
params.newton_solver = 'default';
params.AL_penalty_init = 10;
[out, model] = drspf(model, params);