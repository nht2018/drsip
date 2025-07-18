%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-29 10:53:49
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-01 20:27:18
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-29 10:53:49
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-29 10:55:51
%
%  Copyright (c) 2023, Hantao Nie, Peking University.
%
%%*******************************************************************
%% detect_ublk: this function is modified from detect_ublk.m in SDPT3 
%  by Kim-Chuan Toh, Michael J. Todd, Reha H. Tutuncu
%%*****************************************************************

function [model2, ublkinfo] =  cone_l2u(model)

model2 = model;
K = model.K; At = model.At; c = model.c;
K2 = K; At2 = At; c2 = c;


ublkinfo = cell(length(K), 3);
tol = 1e-14;
%%
numblknew = length(K);
%%
for p = 1:length(K)
    m = size(At{p}, 2);
    cone = K{p};
    if strcmp(cone.type, 'l') || endsWith(cone.type, 'u')
        r = randn(1, m);
        Ap = At{p}'; Cp = c{p};
        ApTr = (r * Ap)';
        [dummy, II] = intersect(ApTr, -ApTr);
        
        if ~isempty(II)
            [tmp, perm] = sort(abs(ApTr(II)));
            idx0 = find(diff(tmp) < tol);
            i1 = II(perm(idx0));
            i2 = II(perm(idx0 + 1));
            n = cone.size;
            Api1 = Ap(:, i1);
            Api2 = Ap(:, i2);
            Cpi1 = Cp(i1)';
            Cpi2 = Cp(i2)';
            
            if ~isempty(i1)
                idxzr = find(abs(Cpi1 + Cpi2) < tol & sum(abs(Api1 + Api2), 1) < tol);
            else
                idxzr = [];
            end
            
            if ~isempty(idxzr)
                i1 = i1(idxzr');
                i2 = i2(idxzr');
                K2{p} = BasicCone('u', length(i1));
                At2{p} = Ap(:, i1)';
                c2{p} = Cp(i1);
                
                fprintf('\n %1.0d linear variables from unrestricted variable.\n', ...
                    2 * length(i1));
                
                
                % the remaining variables is put into a new block
                i3 = setdiff([1:n], union(i1, i2));
                
                if ~isempty(i3)  
                    numblknew = numblknew + 1;
                    K2{numblknew} = BasicCone('l', length(i3));
                    At2{numblknew} = Ap(:, i3)';
                    c2{numblknew} = Cp(i3);
                end
                
                ublkinfo{p, 1} = i1; ublkinfo{p, 2} = i2; ublkinfo{p, 3} = i3;
            end
            
        end
    else
        % do nothing
    end
    
end

model2.K = K2;
model2.At = reshape(At2, [], 1);
model2.c = reshape(c2, [], 1);

end

%%*******************************************************************
