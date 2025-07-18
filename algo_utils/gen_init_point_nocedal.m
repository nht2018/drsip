%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-01-09 17:34:51
%  Description:
%
%  Copyright (c) 2023, Hantao Nie, Peking University.
%

% new update on 2023-11-12
% mu0 is a cell whose size matches that of K
function [mu0, X0, y0, S0, N, n_cones] = gen_init_point1(model, algo, params)


K  = model.K;
X0 = MatCell(length(K));
y0 = zeros(length(model.b), 1);
S0 = MatCell(length(K));
if isfield(params, 'vector_mu') && params.vector_mu
    mu0 = MatCell(length(K));
    mu0_min = 1e-2;
end

Adaggerb = full(algo.Adagger(model.b_int));
AdaggerAc = algo.Adagger(AXfun(model.K, model.At_int, model.c_int));

N = 0;
n_cones = 0;

for p = 1:length(K)
    cone = K{p};
    if strcmp(cone.type, 'l')
        x = Adaggerb{p};
        s = model.c{p} - AdaggerAc{p};
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
            mu0{p} = max(x' * s / sum(cone.size), mu0_min);
        end
        N = N + sum(cone.size);
        n_cones = n_cones + sum(cone.size);
    elseif  strcmp(cone.type, 'b2l')
        len = cone.size;
        x = Adaggerb{p};
        s = [model.c{p}; zeros(len / 2, 1)] - AdaggerAc{p};

        %% shift the linear parts
        min_x = full(min(x));
        min_s = full(min(s));
        delta_x = max(0, -1.5 * min_x);
        delta_s = max(0, -1.5 * min_s);
        x = x + delta_x;
        s = s + delta_s;
        coeff = (x(1: len/2) + x(len/2+1: len)) ./ (cone.params(:, 2) - cone.params(:, 1));
        x(1: len/2) = x(1: len/2) ./ coeff;
        x(len/2+1: len) = x(len/2+1: len) ./ coeff;
        X0{p} = x;
        S0{p} = s;
        if params.vector_mu
            mu0{p} = max(x' * s / sum(cone.size), mu0_min);
        end
        N = N + sum(cone.size);
        n_cones = n_cones + sum(cone.size);
    elseif strcmp(cone.type, 's')
        x = Adaggerb{p};
        s = model.c{p} - AdaggerAc{p};
        
        % shift x and s to the interior of the cone
        [V, D] = eig(x);
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
            mu0{p} = max(inner_product(x, s) / sum(cone.size), mu0_min); % to be modified in the future
        end
        N = N + sum(cone.size);
        n_cones = n_cones + numel(cone.size);
    elseif strcmp(cone.type, 'q')
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
            mu0{p} = max(inner_product(x, s) / sum(cone.size), mu0_min); % to be modified in the future
        end
        N = N + sum(cone.size);
        n_cones = n_cones + numel(cone.size);
    elseif strcmp(cone.type, 'u')
        X0{p} = Adaggerb{p};
        S0{p} = model.c{p} - AdaggerAc{p};
        if params.vector_mu
            mu0{p} = max(sum(X0{p} .* S0{p}) / sum(cone.size), mu0_min);
        end
        N = N + sum(cone.size);
        n_cones = n_cones + sum(cone.size);
    end
end

if ~params.vector_mu
    mu0 = max(inner_product(X0, S0) / N, 1e-2) ;
end


end