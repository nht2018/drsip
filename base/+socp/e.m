%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-11-27 15:07:37
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-27 20:22:04
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = e(cone_size)
    % identity
    ind_head = cumsum(cone_size) - cone_size + 1;
    out = zeros(sum(cone_size), 1);
    out(ind_head) = 1;
end
    