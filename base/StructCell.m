%%  
% we design StructCell to store the Jacobian matrix since for different type of cones, their Jacobian have different form and often have some structure. Hence using a struct to store the Jacobian matrix is more efficient than directly using matrix.
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-07 12:08:11
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-07 12:13:49
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function obj = StructCell(n)
    % StructCell(n) return an object of length n, where n is integer, each block is empty;

    obj = cell(n, 1);
    for i = 1:n
        obj{i} = struct;
    end


end