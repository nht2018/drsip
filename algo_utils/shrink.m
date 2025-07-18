%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-31 17:05:44
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-01 16:30:33
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [X, Z, ind_shrink, model, flag] = shrink(X, Z, model, algo, params)
    % Shrink the quadratic cone to reduce the number of variables
    tol = 1e-7;
    flag = 0;
    ind_shrink = MatCell(length(model.K));
    for p = 1: length(model.K)
        cone = model.K{p};
        ind_shrink{p} = [];
        if strcmp(cone.type, 'q')
            x = X{p};
            ind_head = cumsum(cone.size) - cone.size + 1; 
            ind_zero = setdiff(find(abs(x) < tol) , ind_head); 
            if length(ind_zero)  > 50
                flag = flag + length(ind_zero);
                [count, ~] = histc(ind_zero, [ind_head, length(x) + 1]);
                count = count(1: end - 1);
                model.K{p}.size = cone.size - reshape(count, 1, []);
                X{p}(ind_zero) = [];
                Z{p}(ind_zero) = [];
                % [model, algo] = factAAt_downdate(model.At{p}(ind_zero, :)', model, algo, params);
                % if params.precond_A
                %     precond = algo.AAt;
                %     precond = rmfield(precond, 'time_fact');
                %     model.precond = precond;
                % end
                % model.AAt = model.AAt - model.At{p}(ind_zero, :)' * model.At{p}(ind_zero, :);
                model.At{p}(ind_zero, :) = [];
                model.c{p}(ind_zero) = [];
                ind_shrink{p} = ind_zero;
                fprintf('Shrink %d variables in block %d\n', length(ind_zero), p);
            end
        end
   end

end