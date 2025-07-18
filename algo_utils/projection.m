%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-28 12:05:56
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-09 19:37:28
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-21 19:16:55
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [X, DX] = projection(Z, K)
    % projection mapping onto the cone
    % output
    % X : projection of Z onto the cone
    % DX : derivative of projection mapping


    X = MatCell(length(Z));
    DX = StructCell(length(Z));
    for p =1:length(K)
        cone = K{p};
        if strcmp(cone.type,'l') || endsWith(cone.type, 'l')
            X{p} = max(Z{p}, 0);
            if nargout > 1
                DX{p}.coeff = 0;
                DX{p}.shift = (Z{p} >= 0);
            end
        elseif strcmp(cone.type,'s')
            [V, D] = eig(Z{p});
            d = diag(D);
            X{p} = V * spdiag( max(d, 0) ) * V';
            if nargout > 1
                error('not implemented for sdp cone');
            end
        elseif strcmp(cone.type,'q') || endsWith(cone.type, 'q')
            if nargout == 1
                X{p} = mexprojection_cone_q(Z{p}, cone.size);
            else
                [X{p}, dd, Dsch1, Dsch2, P1, P2, shift] = mexprojection_cone_q(Z{p}, cone.size);
                DX{p}.shift = shift;
                DX{p}.coeff = [Dsch1; Dsch2];
                DX{p}.lr = [P1, P2];
            end
        elseif strcmp(cone.type,'u') 
            X{p} = Z{p};
            if nargout > 1
                DX{p}.coeff = 0;
                DX{p}.shift = ones(size(Z{p}));
            end
        end
    end

end