%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 20:54:19
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-29 17:12:26
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function C = minus(A, B)
    C = elementwiseOperation(A, B, 'minus');
end