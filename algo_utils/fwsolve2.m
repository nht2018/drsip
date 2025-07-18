%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-24 15:19:03
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-10 15:52:56
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [q, info] = fwsolve2(r, info)
    if strcmp(info.solver, 'isidentity') 
        q = r;
    elseif strcmp(info.solver,'isdiag')
        q = r ./ info.sqrtD;
            
    elseif strcmp(info.solver,'ldlchol')
        if ~ (isfield(info, 'L') && isfield(info, 'D') )
            [info.L, info.D] = ldlsplit(info.LD);
            info.D = diag(full(info.D));
        end
        if ~ isfield(info, 'sqrtD')
            info.sqrtD = sqrt(info.D);
        end
        q = (info.L \ r(info.ordering)) ./ info.sqrtD;
    elseif strcmp(info.solver,'lchol')
        q = info.R' \ r(info.ordering);
    elseif strcmp(info.solver,'chol')
        q = info.R' \ r;
    elseif strcmp(info.solver,'spchol')
        q = info.R' \ r(info.ordering);
    else
        error('Unknown solver %s for fwsolve', info.solver);
    end

end
    