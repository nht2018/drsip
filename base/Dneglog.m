%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 16:57:37
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-27 20:17:22
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
%% Derivative of the negative log barrier function
% that is, Dneglog(X, K) = inv(X) if K is not unbounded and 0 otherwise

function result = Dneglog(X, K)
    assert(length(X) == length(K), 'MatCellinv: length(X) ~= length(K)')
    result = MatCell(length(K));
    for p = 1: length(K)
        cone = K{p};
        if strcmp(cone.type, 's')
            result{p} = inv(X{p});
        elseif strcmp(cone.type, 'q')
            result{p} = socp.inv(X{p}, cone.size);
        elseif strcmp(cone.type, 'l')
            result{p} = 1 ./ X{p};
        elseif strcmp(cone.type, 'u')
            result{p} = zeros(size(X{p}));
        end
    end

end
