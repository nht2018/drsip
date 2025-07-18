%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 20:54:42
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-29 17:13:11
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 20:54:42
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-29 17:13:07
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function C = times(A, B)
    C = elementwiseOperation(A, B, 'times');
end