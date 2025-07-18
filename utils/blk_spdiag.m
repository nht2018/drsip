%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-11 21:20:50
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-29 10:23:25
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = blk_spdiag(X, block_size)

%% build a sparse block diagongal matrix with given vector

% input arguments
%     X: a double vector with size [n, 1] or [1, n]
%     block_size: a integer vector with length k, i.e., block_size = [q_1 ,..., q_k]
%         where q_1 + ... + q_k == n

% output
%     a sparse matrix with size = [n, k] if size(X) == [n, 1]
%                              or [k, n] if size(X) == [1, n]

% example
%     >> X = [1: 10]'; block_size = [1 2 3 4];
%     >> blk_spdiag(X, block_size)
%     ans =

%         1     0     0     0
%         0     2     0     0
%         0     3     0     0
%         0     0     4     0
%         0     0     5     0
%         0     0     6     0
%         0     0     0     7
%         0     0     0     8
%         0     0     0     9
%         0     0     0    10

%     >> blk_spdiag(X', block_size)
%     ans =

%         1     0     0     0     0     0     0     0     0     0
%         0     2     3     0     0     0     0     0     0     0
%         0     0     0     4     5     6     0     0     0     0
%         0     0     0     0     0     0     7     8     9    10

    out = mexblk_spdiag(X, block_size);
end