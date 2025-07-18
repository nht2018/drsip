%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-12 11:51:35
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-12 15:03:42
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   


function lhs = linsysinit(At, system_opt)
    m = size(At{1}, 2) ;
    %% initialize lhs 
    lhs = struct();

    lhs.mat11 = sparse(m, m);
    if strcmp(system_opt, 'sparse') || strcmp(system_opt, 'sparse_psd')
        lhs.mat21_den = sparse(0, m);
        lhs.mat21_u = sparse(0, m);
    elseif strcmp(system_opt, 'dense') || strcmp(system_opt, 'augmented')
        lhs.mat21 = sparse(0, m);
    else
        error("unknown system_opt: %s", system_opt);
    end
    lhs.mat31 = sparse(0, m);
    lhs.mat32 = [];
    lhs.mat22_diag = [];
    lhs.mat33_diag = [];

end