%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-18 17:07:34
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-11 21:04:17
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function Y = AAtfun(At)
% AAtfun  Compute A*A' for a given matrix A
if iscell(At) 
    Y = sparse(size(At{1}, 2), size(At{1}, 2));
    for p = 1: length(At)
        Y = Y + At{p}' * At{p};
    end
else % At is a single matrix
    Y = At' * At;
end