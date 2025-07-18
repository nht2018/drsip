%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-11-12 16:11:04
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-24 20:00:01
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = min_all(obj)
    % compute the sum of all the blocks
    if length(obj) == 0
        out = 0;
        return
    end

    out = min(cellfun(@(x) min(x(:)), obj));

end