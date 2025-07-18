%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-23 17:19:28
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-01-09 16:18:54
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-11 21:04:02
%  Copyright (c) 2023, Hantao Nie, Peking University.
%

%% comput Aty = A'*y

function Aty = Atyfun(K, At, y)

Aty = MatCell(length(K));
for p = 1:length(K)
    cone = K{p};
    n = sum(cone.size);
    if strcmp(cone.type,'s')
        if (isempty(At{p}))
            Aty{p} = sparse(n,n);
        else
            Aty{p} = mexsmat({cone.type, cone.size}, At{p}*y);
        end
    else
        if (isempty(At{p}))
            Aty{p} = zeros(n,1);
        else
            Aty{p} = At{p}*y;
        end
    end
end

end
