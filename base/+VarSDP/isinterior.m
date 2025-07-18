%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-02 17:55:21
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-22 15:39:59
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
% checks if X is in the interior of the cone K

function [flag] = isinterior(K, X)
    flag = true;
    for p = 1: length(K)
        cone = K{p};
        if strcmp(cone.type, 'l') || endsWith(cone.type, 'l')
            if min(X{p}) < 0
                flag = false;
                return
            end
        elseif strcmp(cone.type, 'q') || endsWith(cone.type, 'q')
            if min(socp.det(X{p}, cone.size)) < 0
                flag = false;
                return
            end
        elseif strcmp(cone.type, 's')
            if eigs(X{p}, 1, 'smallestabs') < 0
                flag = false;
                return
            end
        elseif strcmp(cone.type, 'u')
            % pass
        else
            error("unknown cone type")
        end
    end
end