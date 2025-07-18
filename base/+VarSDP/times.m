%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-02 09:49:38
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-02 09:55:02
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = times(K, X, Y)
    out = MatCell(length(K));
    for p = 1: length(K)
        cone = K{p};
        if strcmp(cone.type, 'q')
            out{p} = socp.times(X{p}, Y{p}, cone.size) ;
        elseif strcmp(cone.type, 'l') || strcmp(cone.type, 'u')
            out{p} = X{p} .* Y{p} ;
        else
            error("not implemented yet")
        end
    end
end