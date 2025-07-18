%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-11-27 20:19:49
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-27 20:23:04
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = sqrtdet(X, cone_size)
    out = mexcone_q_ops(X, 'sqrtdet', cone_size);
end