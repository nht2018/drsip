%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-23 17:19:28
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-01-04 15:26:37
%
%  Copyright (c) 2024, Hantao Nie, Peking University.
%
%% compute AX = A * X

%% A * X
function AX = AXfun(K, At, X)
m = size(At{1}, 2);
AX = zeros(m, 1);
for p = 1:length(K)
    cone = K{p};
    if isempty(At{p})
        continue;
    end
    Xp = X{p};
    if strcmp(cone.type, 's')
        % if isa(X{p}, "Var_sdp")  % X{p} is Var_sdp object
        %    x = mysvec(X{p});
        %    AX = AX + (x'*At{p})';
        % else % X{p} is one single matrix or vector.
        if length(cone.size) > 1
            Xp = sparse(Xp);
        end
        
        if size(X{p}, 1) == sum(cone.size .* (cone.size + 1) / 2) % X{p} is in vector form
            AX = AX + At{p}' * Xp;
        elseif size(X{p}, 1) == sum(cone.size) % X{p} is in matrix form
            x = mysvec(cone, Xp);
            AX = AX + At{p}' * x;
        else
            error('AXfun: dimension mismatch for block %i', p);
        end
    else
        AX = AX + At{p}' * Xp;
    end
end
end

