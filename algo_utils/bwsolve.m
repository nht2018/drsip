%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-24 15:30:50
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-06 22:52:04
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function q = bwsolve(factor,r)
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
            q(factor.ordering,1) = factor.L' \ (factor.sqrtD \ r);
        elseif strcmp(factor.solver,'lchol')
            q(factor.ordering,1) = factor.R \ r;
        elseif strcmp(factor.solver,'chol')
            % q(factor.ordering,1) = mextriang(factor.R, r, 1);
            q = factor.R \ r;
        elseif strcmp(factor.solver,'spchol')
            % q(factor.ordering,1) = mexbwsolve(factor.Rt, r);
            q(factor.ordering,1) = factor.R \ r;
        end
    end
end