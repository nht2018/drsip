%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-07 11:48:08
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-09 19:36:48
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
%% Now we support only affine and affine_inv 
% Given a ProjJac object jac. jac can be seen as a mappping.
% affine means 
%   out = scale * jac  + shift * identity
% affine_inv means
%   out = (scale * jac  + shift * identity)^{-1}
% out is also a ProjJac object

function out = Jac_ops(jac, K, operand, scale, shift)
    assert(strcmp(operand, 'affine') || strcmp(operand, 'affine_inv'))
    
    % check jac and K are compatible
    % to be implemented

    out = StructCell(length(K));
    


    % transformation on coeff, Dsch2 and shift
    if strcmp(operand, 'affine')
        for p = 1: length(K)
            cone = K{p};
            if strcmp(cone.type, 's') 
                out{p}.coeff = scale * jac{p}.coeff;
                out{p}.shift = scale * jac{p}.shift + shift;
                out{p}.eigs = jac{p}.eigs;
            elseif  strcmp(cone.type, 'q') || endsWith(cone.type, 'q')
                out{p}.coeff = scale * jac{p}.coeff;
                out{p}.shift = scale * jac{p}.shift + shift;
                out{p}.lr = jac{p}.lr;
            elseif strcmp(cone.type, 'l') || endsWith(cone.type, 'l') ||strcmp(cone.type, 'u') 
                out{p}.shift = scale * jac{p}.shift + shift;
            end
        end

        % % check
        % r = MatCell(length(K));
        % for p = 1: length(K)
        %     cone = K{p};
        %     r{p} = rand(sum(cone.size), 1);
        % end
        % temp1 = Jac_lmut(K, jac, r);
        % temp1 = scale * temp1 + shift * r;
        % temp2 = Jac_lmut(K, out, r);
        % fprintf('check affine: %e\n', norm(temp1 - temp2) / norm(r));

    elseif strcmp(operand, 'affine_inv')
        for p = 1: length(K)
            cone = K{p};
            if strcmp(cone.type, 's')
                out{p}.shift = 1 ./ (scale * jac{p}.shift  + shift);
                out{p}.shift(abs(out{p}.shift) > 1e16) = 1e16 * sign(out{p}.shift(abs(out{p}.shift) > 1e16));
                out{p}.coeff = 1 ./ (scale * ( jac{p}.coeff + jac{p}.shift ) + shift) - out{p}.shift;
                out{p}.eigs = jac{p}.eigs;
            elseif strcmp(cone.type, 'q') || endsWith(cone.type, 'q') 
                assert(size(jac{p}.shift, 1) * 2  == size(jac{p}.coeff, 1));  
                out{p}.shift = 1 ./ (scale * jac{p}.shift  + shift); % size = [length(cone.size), 1]
                out{p}.shift(abs(out{p}.shift) > 1e16) = 1e16 * sign(out{p}.shift(abs(out{p}.shift) > 1e16));
                out{p}.coeff = - [out{p}.shift; out{p}.shift] .* (scale * jac{p}.coeff ./ (scale * jac{p}.coeff + scale * [jac{p}.shift; jac{p}.shift] + shift) ); % size = [length(cone.size) * 2, 1]
                out{p}.coeff(abs(jac{p}.coeff) < 1e-16 ) = 0;
                out{p}.coeff(abs(out{p}.coeff) > 1e16) = 1e16 * sign(out{p}.coeff(abs(out{p}.coeff) > 1e16));
                out{p}.lr = jac{p}.lr;
            elseif strcmp(cone.type, 'l') || endsWith(cone.type, 'l')  || strcmp(cone.type, 'u')
                out{p}.shift = 1 ./ (scale * jac{p}.shift + shift);
                out{p}.shift(abs(out{p}.shift) > 1e16) = 1e16 * sign(out{p}.shift(abs(out{p}.shift) > 1e16));
            end
        end

        % % check
        % r = MatCell(length(K));
        % for p = 1: length(K)
        %     cone = K{p};
        %     r{p} = rand(sum(cone.size), 1);
        % end
        % jac1 = Jac_ops(jac, K, 'affine', scale, shift);
        % temp1 = Jac_lmut(K, jac1, r);
        % temp2 = Jac_lmut(K, out, temp1);
        % % for p = 1: length(K)
        % %     fprintf('check affine_inv block %i, type %s: %e\n', p, K{p}.type, norm(temp2{p} - r{p}) / norm(r{p}));
        % % end
        % fprintf('check affine_inv: %e\n', norm(r - temp2) / norm(r));


    end


    
end