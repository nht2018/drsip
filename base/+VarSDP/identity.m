%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-02 09:49:38
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-02 22:15:00
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = identity(K)
    out = MatCell(length(K));
    for p = 1: length(K)
        cone = K{p};
        if strcmp(cone.type, 'q')
            out{p} = socp.e(cone.size);
        elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'u')
            out{p} = ones(sum(cone.size), 1);
        elseif strcmp(cone.type, 's')
            out{p} = speye(sum(cone.size), sum(cone.size));
        end
    end
end