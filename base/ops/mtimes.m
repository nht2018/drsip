%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 21:04:46
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-29 17:24:30
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   

function C = mtimes(A, B)  
    C = elementwiseOperation(A, B, 'mtimes');
end