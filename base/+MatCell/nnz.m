%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 21:02:48
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-24 15:30:46
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function out = nnz(obj)
    % sum of non-zero elements in cell array
    out = sum(cellfun(@nnz, obj));
end