%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-01-29 20:30:56
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-22 15:35:05
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
%% 
% std_model has the form as
% min <c, x> + c0
% s.t. A * x = b
%      x \in K
% 
% the model struct has fields as
% At: MatCell
% c: MatCell
% b: double vector
% K: Cone
% c0: double, optional


function flag = check_std_model(model)
    % check fields
    if ~all(isfield(model, {'At', 'c', 'b', 'K'}))
        flag = false;
        return;
    end

    if length(model.K) ~= length(model.At) || length(model.K) ~= length(model.c)
        flag = false;
        return;
    end

    % check dimensions
    for p = 1:length(model.K)
        % check number of variables
        cone = model.K{p};
        if strcmp(cone.type, 's')
            if size(model.At{p}, 1) ~= sum(cone.size .* (cone.size + 1) / 2) || size(model.c{p}, 1) ~= sum(cone.size) || size(model.c{p}, 2) ~= sum(cone.size)
                flag = false;
                return;
            end

        else
            if size(model.At{p}, 1) ~= length(model.c{p})
                flag = false;
                return;
            end
        end

        % check number of constraints
        if size(model.At{p}, 2) ~= length(model.b)
            flag = false;
            return;
        end
    end



    flag = true;
end