%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-02 09:49:38
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-02 09:56:18
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = inv(K, X)
    out = MatCell(length(K));
    for p = 1: length(K)
        cone = K{p};
        if strcmp(cone.type, 'q')
            out{p} = socp.inv(X{p}, cone.size) ;
        elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'u') || strcmp(cone.type, 'b2l')
            out{p} = 1 ./ X{p};
        else
            error("not implemented yet")
        end
    end
end