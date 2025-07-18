%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-06 13:35:51
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   

% new update on 2023-11-12
% mu0 is a cell whose size matches that of K
function [Z0, mu0, X0, S0, n_cones] = gen_init_point1(model, algo, params)
    K  = model.K;
    X0 = MatCell(length(K));
    S0 = MatCell(length(K));
    if params.vector_mu
        mu0 = MatCell(length(K));
        mu0_min = 1e-2;
    end
    
    Adaggerb = full(algo.Adagger(model.b_int));
    AdaggerAc = algo.Adagger(algo.AXmap_int(model.c_int));


    n_cones = 0;
    for p = 1:length(K)
        cone = K{p};
        if strcmp(cone.type, 'l') || endsWith(cone.type, 'l')
            if strcmp(cone.type, 'b2l')
                x = Adaggerb{p};
                s = [model.c{p}; zeros(cone.size / 2, 1)] - AdaggerAc{p};
            else
                x = Adaggerb{p};
                s = model.c{p} - AdaggerAc{p};
            end
            %% shift the linear parts
            min_x = full(min(x));
            min_s = full(min(s));
            delta_x = max(0, -1.5 * min_x);
            delta_s = max(0, -1.5 * min_s);
            x = x + delta_x;
            s = s + delta_s;
            X0{p} = x;
            S0{p} = s;
            if params.vector_mu
                mu0{p} = max(x' * s / algo.AL_penalty / length(x), mu0_min);
            end
            n_cones = n_cones + length(x);
        elseif strcmp(cone.type, 's')
            x = Adaggerb{p};
            s = model.c{p} - AdaggerAc{p};
        
            % shift x and s to the interior of the cone
            [V, D] = eig(full(x));
            d = diag(D);
            min_d = min(d);
            delta_d = max(0, -1.5 * min_d);
            x = V * spdiag(d + delta_d) * V';
            [V, D] = eig(full(s));
            d = diag(D);
            min_d = min(d);
            delta_d = max(0, -1.5 * min_d);
            s = V * spdiag(d + delta_d) * V';
            X0{p} = x;
            S0{p} = s;
            if params.vector_mu
                mu0{p} = max(inner_product(x, s) / algo.AL_penalty / sum(cone.size), mu0_min); % to be modified in the future
            end
            n_cones = n_cones + numel(cone.size);
        elseif strcmp(cone.type, 'q') || endsWith(cone.type, 'q')
            x = Adaggerb{p};
            s = model.c{p} - AdaggerAc{p};
            % shift x and s to the interior of the cone
            cone = model.K{p};
            ind_head = cumsum(cone.size) - cone.size + 1; 
            x(ind_head) = max(x(ind_head), 1.5 * sqrt(soc_ops(full(x), 'norm', cone.size) .^ 2 - x(ind_head) .^ 2) );
            s(ind_head) = max(s(ind_head), 1.5 * sqrt(soc_ops(full(s), 'norm', cone.size) .^ 2 - s(ind_head) .^ 2) );
            X0{p} = x;
            S0{p} = s;
            if params.vector_mu
                mu0{p} = max(inner_product(x, s) / algo.AL_penalty / sum(cone.size), mu0_min); % to be modified in the future
            end
            n_cones = n_cones + numel(cone.size);
        elseif strcmp(cone.type, 'u')
            X0{p} = Adaggerb{p};
            S0{p} = model.c{p} - AdaggerAc{p};
            if params.vector_mu
                mu0{p} = max(sum(X0{p} .* S0{p}) / algo.AL_penalty / sum(cone.size), mu0_min);
            end
            n_cones = n_cones + sum(cone.size);
        end
    end

    Z0 = X0 - S0;
    if ~params.vector_mu
        mu0 = max(inner_product(X0, S0)/ algo.AL_penalty / n_cones, 1e-2) ;
    end


end