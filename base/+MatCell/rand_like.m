%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 21:06:57
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-24 12:57:51
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function Y = rand_like(X)
    Y = cell(size(X));
    for i = 1:length(X)
        Y{i} = rand(size(X{i}));
    end
end