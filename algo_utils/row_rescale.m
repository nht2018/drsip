%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-24 16:59:40
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-24 19:35:52
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [model, scale_A] = row_rescale(model)

At = model.At;
b = model.b;
m = size(model.b, 1);
row_normA = zeros(m, 1);
for p = 1: length(model.K)
    row_normA = row_normA + (vecnorm(At{p}, 2, 1)') .^ 2;
end
row_normA = sqrt(row_normA);

% empty_row = find(row_normA == 0);
% if ~isempty(empty_row)
%     fprintf("remove %d empty rows in A.\n", length(empty_row) );
%     for p = 1: length(model.K)
%         At{p}(:, empty_row) = [];
%     end
%     b(empty_row) = [];
%     row_normA(empty_row) = [];
% end
row_normA(row_normA == 0) = 1;


for p = 1: length(model.K)
    At{p} = At{p} * spdiag(1 ./ row_normA);
end

b = b ./ row_normA;
scale_A = row_normA;

model.At = At;
model.b = b;

end


