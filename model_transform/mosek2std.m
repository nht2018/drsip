%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-12 12:06:29
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-29 21:31:54
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function model = mosek2std(mosekmodel)
    %% transform a mosek model to standard model

    % input problem is
    % min c' * x + cfix
    % s.t.  blc <= a * x <= buc
    %     f * x + g in accs
    %     blx <= x <= bux

    model = mosek2gen(mosekmodel);

    % gen2std my model to be standard form
    model = gen2std(model);

 

    
   
end