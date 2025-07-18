%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-06 16:22:06
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-06 16:25:23
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function print_struct(fid, s, indent)
    if nargin < 3
        indent = '';
    end
    fn = fieldnames(s);
    for i = 1:numel(fn)
        if isstruct(s.(fn{i}))
            fprintf(fid, indent);
            fprintf(fid, '%s:\n', fn{i});
            print_struct(fid, s.(fn{i}), [indent '\t']);
        else
            fprintf(fid, indent);
            fprintf(fid, '%s: %s\n', fn{i}, num2str(s.(fn{i})));
        end
    end
end