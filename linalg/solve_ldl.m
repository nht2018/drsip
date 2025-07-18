%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-03-05 15:52:40
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-05 15:52:53
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [d] = solve_ldl(rhs, info)
    assert(isfield(info, 'L'));
    assert(isfield(info, 'D'));
    d = info.L \ rhs; d = info.D \ d; d = info.L' \ d;
end