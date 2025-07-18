%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-29 11:09:15
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-06 15:00:57
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [model] = cone_u2l(model)
for p = 1: length(model.K)
    cone = model.K{p};
    if strcmp(cone.type, 'u')
        model.K{p} = BasicCone('u2l', 2 * sum(cone.size));
        model.At{p} = [model.At{p}; -model.At{p}];
        model.c{p} = [model.c{p}; - model.c{p}];
    end
end

end