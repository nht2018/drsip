%%  
%  compute the Arrow shape matrix of a vector Arw(X) or compute the matrix product Arw(X) * Y
%  usage:
%     out = Arw(X, 'mat', cone_size)
%     out = Arw(X, Y, cone_size)

%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-11-27 12:41:01
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-27 20:57:36
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = Arw(X, Y, cone_size)

    if strcmp(Y, 'mat')
        J = socp.J(cone_size);
        e = socp.e(cone_size);
        ind_head = cumsum(cone_size) - cone_size + 1;
        out = - spdiag(repelem(X(ind_head), cone_size) .* J ) + blk_spdiag(e, cone_size) * blk_spdiag(X, cone_size)' +  blk_spdiag(X, cone_size) * blk_spdiag(e, cone_size)';
    else
        out = socp.times(X, Y, cone_size);
    end
end