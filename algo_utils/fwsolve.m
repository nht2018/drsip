%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-24 15:19:03
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-10 15:52:56
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function q = fwsolve(factor, r)
    if (isfield(factor, 'isidentity') && factor.isidentity)
        q = r;
    else
        if strcmp(factor.solver,'ldlchol')
            if ~ (isfield(factor, 'L') && isfield(factor, 'D') )
                [factor.L, factor.D] = ldlsplit(factor.LD);
            end
            if ~ isfield(factor, 'sqrtD')
                factor.sqrtD = sqrt(factor.D);
            end
            q = factor.sqrtD \ (factor.L \ r(factor.ordering));
        elseif strcmp(factor.solver,'lchol')
            q = factor.Rt \ r(factor.ordering);
        elseif strcmp(factor.solver,'chol')
            % q = mextriang(factor.R, r(factor.ordering), 2);
            q = factor.Rt \ r;
        elseif strcmp(factor.solver,'spchol')
            % q = mexfwsolve(factor.R, r(factor.ordering,1));
            q = factor.Rt \ r(factor.ordering);
        end
    end
end
    