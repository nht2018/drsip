%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-02-26 20:46:59
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-26 20:54:16
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [d] = solve_chol(rhs, info)
    assert(isfield(info, 'R'));
    d = info.R \ ( info.R' \ rhs );
end