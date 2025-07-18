%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 20:56:54
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-29 17:12:59
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function C = sparse(A)
    C = unaryOperation(A, 'sparse');
end