%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-18 10:54:06
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-08 20:55:50
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function sigPow = sigPow_update(normF)
    if normF < 1e-2
        sigPow      = 0.4;
    else
        sigPow      = 0.5;
    end
end
