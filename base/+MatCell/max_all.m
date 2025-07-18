%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-11-12 16:11:04
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-24 20:00:05
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = max_all(obj)
    % compute the sum of all the blocks
    if length(obj) == 0
        out = 0;
        return
    end

    out = max(cellfun(@(x) max(x(:)), obj));

end