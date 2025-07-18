%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-05 22:22:46
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [Z0, mu0, X0, S0, n_cones] = gen_init_point2(model, algo, params)

%% This heuristic is motivated by SDPT3


    K  = model.K;
    X0 = MatCell(length(K));
    S0 = MatCell(length(K));
    At = model.At_int;
    b  = model.b_int;
    c = model.c_int;
    % purturb = 1e-10; % from SDPT3
    purturb = 0;

    n_cones = 0;
    for p = 1:length(K)
        cone = K{p};
        if strcmp(cone.type, 'l') || endsWith(cone.type, 'l')
            len = length(c{p});
            normc = 1 + norm(c{p});
            normA = 1 + sqrt(sum(At{p}.^2, 1)'); %norm of each row of A
            const = 10;
            constX = max([const, sqrt(len), max(sqrt(len) * (1 + abs(b))./ normA)]); 
            constS = max([const,sqrt(len), max(normA), normc]);
            x = constX*(1+purturb*rand(len, 1));
            s = constS*(1+purturb*rand(len, 1));

            X0{p} = x;
            S0{p} = s;
            n_cones = n_cones + sum(cone.size);
        elseif strcmp(cone.type, 's')
            normAni = [];
            n = sum(cone.size);
            x = sparse(n,n); 
            s = sparse(n,n);
            ss = [0, cumsum(cone.size)];
            tt = [0, cumsum(cone.size.*(cone.size+1)/2)];
            for i = 1:length(cone.size)
               if ~isempty(At{p})
                  pos = [tt(i)+1 : tt(i+1)];
                  Ai = At{p}(pos,:);
                  normAni = 1+sqrt(sum(Ai.*Ai));
               end
               pos = [ss(i)+1 : ss(i+1)];  
               ni = length(pos);
               tmp = c{p}(pos,pos);
               normCni = 1+sqrt(sum(sum(tmp.*tmp)));
               const  = 10; %%--- old: const = 1; 
               constX = max([const,sqrt(ni),ni*( (1 + abs(b))./ normAni)]); 
               constZ = max([const,sqrt(ni),normAni,normCni]);
               x(pos,pos) = constX*spdiags(1+purturb * rand_psd(ni),0,ni,ni);
               s(pos,pos) = constZ*spdiags(1+purturb * rand_psd(ni),0,ni,ni);
            end

            X0{p} = x;
            S0{p} = s;
            n_cones = n_cones + numel(cone.size);
        elseif strcmp(cone.type, 'q') || endsWith(cone.type, 'q')
            len = sum(cone.size);
            ind_head = cumsum(cone.size) - cone.size + 1;
            n_subblk = length(cone.size);
            normc = 1+norm(c{p});
            normA = 1+sqrt(sum(At{p} .^ 2, 1)' ); % norm of each row of A
            x = zeros(len,1);
            s = zeros(len,1);
            x(ind_head) = max(1, max((1 + abs(b))./normA) )*sqrt(cone.size') ;
            s(ind_head) = max([sqrt(cone.size); max(max(normA) , normc)*ones(1, n_subblk)])';
            x(ind_head) = x(ind_head).*(1+purturb*rand(n_subblk, 1));
            s(ind_head) = s(ind_head).*(1+purturb*rand(n_subblk, 1));

            X0{p} = x;
            S0{p} = s;
            n_cones = n_cones + numel(cone.size);
        elseif strcmp(cone.type, 'u')
            X0{p} = zeros(sum(cone.size), 1);
            S0{p} = zeros(sum(cone.size), 1);
            n_cones = n_cones + sum(cone.size);
        end
    end

    Z0 = X0 - S0;
    mu0 = max(inner_product(X0, S0) / algo.AL_penalty / n_cones, 1e-2);


end