%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-24 15:30:50
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-06 22:52:04
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [q, info] = bwsolve2(r, info)
    if strcmp(info.solver,'isidentity')
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
        q = zeros(size(r));
        q(info.ordering, :) = info.L' \ (r ./ info.sqrtD);
    elseif strcmp(info.solver,'lchol')
        q = zeros(size(r));
        q(info.ordering, :) = info.R \ r;
    elseif strcmp(info.solver,'chol')
        q = info.R \ r;
    elseif strcmp(info.solver,'spchol')
        q = zeros(size(r));
        q(info.ordering, :) = info.R \ r;
    else
        error('Unknown solver %s for bwsolve', info.solver);
    end
end