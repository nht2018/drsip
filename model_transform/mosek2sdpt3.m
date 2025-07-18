%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-04 20:22:26
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function model = mosek2sdpt3(model)

    model = mosek2std(model);

    model.C = model.c;
    model = rmfield(model, 'c');
    model.blk = Cone.toblk(model.K);
    model = rmfield(model, 'K');
 

    
   
end