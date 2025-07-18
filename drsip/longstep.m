%%  
% compute the long step in neighborhood
% that is, the maximum alpha such that
% X + alpha * dX \in K
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-02 15:36:45
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-01-09 20:09:46
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [stepsize] = longstep(X, dX, model, algo, params)

    K = model.K;
    stepsize = 1;
    for p = 1: length(K)
        cone = K{p};
        if strcmp(cone.type, 'l')   
            idx = find(dX{p} < 0); 
            if ~isempty(idx)
                alpha = min(-X{p}(idx)./dX{p}(idx));  
                stepsize = min([stepsize, alpha]);
            else 
                % pass
            end
        elseif strcmp(cone.type, 'b2l')  
            idx = find(dX{p} < 0); 
            if ~isempty(idx)
                alpha = min(-X{p}(idx)./dX{p}(idx));  
                stepsize = min([stepsize, alpha]);
            else 
                alpha = 1e12;
            end
            % fprintf("step of box: %.4e\n", alpha);
        elseif strcmp(cone.type, 'u')
            % pass
        elseif strcmp(cone.type, 'q')
            a = socp.det(dX{p}, cone.size);
            b = 2 * socp.J_inner_product(X{p}, dX{p}, cone.size);
            c = socp.det(X{p}, cone.size);
            delta = b .^ 2 - 4 * a .* c;
            idx = find(delta > 0 & min(a, b) < 0); 
            steptmp = 1e12 * ones(length(cone.size),1); 
            if ~isempty(idx)
                steptmp(idx) = -(b(idx)+sqrt(delta(idx)))./ (2 * a(idx) );       
            end
            idx = find(abs(a) < 1e-12 & b < 0); 
            if ~isempty(idx)
                steptmp(idx) = -c(idx)./ b(idx); 
            end
            ind_head = cumsum(cone.size) - cone.size + 1;
            dX0 = dX{p}(ind_head); 
            X0 = X{p}(ind_head); 
            idx = find(dX0 < 0 & X0 > 0); 
            if ~isempty(idx)
                steptmp(idx) = min(steptmp(idx), -X0(idx) ./ dX0(idx)); 
            end
            alpha = min(steptmp);
            stepsize = min([stepsize, alpha]);
        end


    end
end