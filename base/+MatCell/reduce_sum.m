%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-23 17:15:00
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-25 22:29:35
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   

function out = reduce_sum(obj)
    % compute the sum of all the blocks
    if length(obj) == 0
        out = [];
        return
    end

    out = obj{1};
    for k = 2: length(obj)
        out = out + obj{k};
    end



end