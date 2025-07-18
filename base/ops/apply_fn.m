%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-12 09:50:29
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-29 17:43:42
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function C = apply_fn(A, B)
    if iscell(A) && iscell(B)
        C = cellfun(@(a, b) a(b), A, B, 'UniformOutput', false);
    elseif iscell(A) 
        C = cellfun(@(a) a(B), A, 'UniformOutput', false);
    elseif iscell(B)
        C = cellfun(@(b) A(b), B, 'UniformOutput', false);
    else
        C = A(B);
    end
end