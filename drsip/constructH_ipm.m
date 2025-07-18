%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-01-14 10:07:42
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-01-29 11:59:10
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [algo] = constructH(X, y, S, mu, model, algo, params)
K = model.K;
H = StructCell(length(K));
maxH = 1e16;
minH = 1e-16;
for p = 1:length(K)
    cone = K{p};
    if strcmp(cone.type,'l') || strcmp(cone.type,'b2l')
        H{p}.shift = X{p} ./ S{p} ;
        H{p}.shift = min(H{p}.shift, maxH);
        H{p}.shift = max(H{p}.shift, minH);
    elseif strcmp(cone.type,'u') 
        H{p}.shift = 0;
    elseif strcmp(cone.type, 's')
        % H.lmut{p} = @(r_) 1 / ( (1 + sigma) * mu) * X{p} * r_ * X{p}; % here r_ is a matrix
        error("not implemented yet");
    elseif strcmp(cone.type,'q') || strcmp(cone.type,'r2q')
        if strcmp(algo.newton.scaling, 'HKM')
            k = length(cone.size);
            x = full(X{p});
            s = full(S{p});
            sinv = socp.inv(s, cone.size);
            dets = socp.det(s, cone.size);
            J = socp.J(cone.size);
            % H{p}.shift = - repelem(blk_sum(x .* s, cone.size, 1, false) ./ dets, cone.size, 1) .* J;
            % H{p}.coeff = [sparse(k, k), speye(k, k);
            %     speye(k, k), sparse(k, k)];
            % H{p}.invcoeff = [sparse(k, k), speye(k, k);
            %     speye(k, k), sparse(k, k)];
            % H{p}.lr = [blk_spdiag(x, cone.size), blk_spdiag(sinv, cone.size)] ;
            
            H{p}.shift = - repelem(blk_sum(x .* s, cone.size, 1, false) ./ dets, cone.size, 1) .* J;
            H{p}.coeff = [sparse(k, k), spdiag(1 ./ dets);
                        spdiag(1 ./ dets), sparse(k, k)];
            H{p}.invcoeff = [sparse(k, k), spdiag(dets);
                        spdiag(dets), sparse(k, k)];
            H{p}.lr = [blk_spdiag(x, cone.size), blk_spdiag(sinv .* repelem(dets, cone.size, 1), cone.size)] ;


            % % check H
            % r1 = rand(sum(cone.size), 1);
            % % r2 = H{p}.shift .* r1 + H{p}.lr * (H{p}.coeff * (H{p}.lr' * r1));
            % r2 = H{p}.shift .* r1 + blk_spdiag(x, cone.size) * blk_spdiag(sinv, cone.size)' * r1 + blk_spdiag(sinv, cone.size) * blk_spdiag(x, cone.size)' * r1;
            % P = socp.sqrt(S{p}, cone.size) ;
            % Q = socp.sqrt(socp.inv(S{p}, cone.size), cone.size) ;
            % tilde_x = socp.Q(P, X{p}, cone.size);
            % r3 = socp.Q(Q, r1, cone.size);
            % r3 = socp.times(tilde_x, r3, cone.size);
            % r3 = socp.Q(Q, r3, cone.size);
            % fprintf("check H res = %e\n", norm(r2 - r3) );
        elseif strcmp(algo.newton.scaling, 'none')
            error("not implemented yet");
            % H = Arw(S) ^{-1} * Arw(X) 
            H{p}.shift = 0;
            H{p}.coeff = 0;
            H{p}.lr = 0;
        else
            error("not implemented yet");
        end
    end
end




%% handle box constraints
H_rd = H;
for p = 1:length(K)
    cone = K{p};
    if strcmp(cone.type, 'b2l')
        % reduce the dim from cone.size to cone.size / 2
        h1 = H{p}.shift(1: cone.size / 2);
        h2 = H{p}.shift(cone.size / 2 + 1: cone.size);
        H_rd{p}.shift = h1 .* h2 ./ (h1 + h2);
    end
end


algo.newton.H = H;
algo.newton.H_rd = H_rd;

maxH = 0;
minH = inf;
for p = 1:length(K)
    if ~ strcmp(K{p}.type, 'u')
        maxH = max(maxH, max(abs(H{p}.shift)));
        minH = min(minH, min(abs(H{p}.shift)));
    end
end
% minH
% maxH
algo.newton.condH = maxH / minH;
end