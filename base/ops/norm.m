%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-11 20:57:13
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-17 22:43:59
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function out = norm(varargin)
    if nargin == 1 && iscell(varargin{1}) 
        out = sqrt(sum(cellfun(@(x) norm(x, 'fro'), varargin{1}) .^ 2));
    else
        out = builtin('norm', varargin{:});
    end
end