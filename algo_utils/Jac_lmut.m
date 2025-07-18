%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-28 12:05:23
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-22 15:48:06
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
% compute V * Z
% where V is Jacobian of proximal operator or projection operator

function out = Jac_lmut(K, V, Z)
    out = MatCell(length(K));
    for p= 1: length(K)
        cone = K{p};
        if strcmp(cone.type, 'l') || endsWith(cone.type, 'l')
            out{p} = V{p}.shift .* full(Z{p});
        elseif strcmp(cone.type, 's')
            if isfield(V{p}, 'lmut')
                out{p} = V{p}.lmut(Z{p});
            else
                out{p} = V{p}.eigs * (V{p}.coeff .* (V{p}.eigs' * Z{p} * V{p}.eigs)) * V{p}.eigs' + V{p}.shift .* Z{p};
            end
        elseif strcmp(cone.type,'q') || endsWith(cone.type, 'q')
            % if size(V{p}.shift, 1) == length(cone.size) then duplicate
            if length(V{p}.shift) == length(cone.size)
                V{p}.shift = repelem(V{p}.shift, cone.size, 1);
            end
            if size(V{p}.coeff, 2 ) == 1
                V{p}.coeff = spdiag(V{p}.coeff);
            end
            if size(V{p}.lr, 2) < size(V{p}.coeff, 1)
                lr_cell = arrayfun(@(j) blk_spdiag(V{p}.lr(:, j), cone.size), 1:size(V{p}.lr, 2), 'UniformOutput', false);
                lr_ = horzcat(lr_cell{:});
                % lr_ = [blk_spdiag(V{p}.lr(:, 1), cone.size), blk_spdiag(V{p}.lr(:, 2), cone.size)];
            else
                lr_ = V{p}.lr;
            end
            out{p} = V{p}.shift .* full(Z{p}) + full(lr_ * (V{p}.coeff * full(lr_' * Z{p}) ));
        elseif strcmp(cone.type,'u') 
            out{p} = V{p}.shift .* full(Z{p});
        end
    end
end