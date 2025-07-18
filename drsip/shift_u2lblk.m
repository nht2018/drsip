%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-29 11:18:26
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-29 11:40:17
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [X, y, S] = shift_cone_u2l(X, y, S, mu, step_d, u2lflag, model, algo, params)
    dres =  algo.Atymap(y) + S - model.c;
    dinf = norm(dres) / (1 + norm(model.c));
    for p = 1: length(model.K)
        if u2lflag(p)
            cone = model.K{p};
            len = cone.size / 2;
            idx1 = 1: len;
            idx2 = len + 1: 2 * len;
            xtmp = min(X{p}(idx1), X{p}(idx2));
            alpha = 0.8;
            X{p}(idx1) = X{p}(idx1) - alpha * xtmp;
            X{p}(idx2) = X{p}(idx2) - alpha * xtmp;
            
            if (mu < 1e-4) % % old: (mu < 1e-7)
                S{p} = 0.5 * mu ./ max(1, X{p}); % % good to keep this step
            else
                stmp = min(1, max(S{p}(idx1), S{p}(idx2)));
                
                if (  dinf > 1e-4 & step_d < 0.2)
                    beta = 0.3;
                else
                    beta = 0.0;
                end
                
                %% important to set beta = 0 at later stage.
                S{p}(idx1) = S{p}(idx1) + beta * stmp;
                S{p}(idx2) = S{p}(idx2) + beta * stmp;
            end
        end
    end
end