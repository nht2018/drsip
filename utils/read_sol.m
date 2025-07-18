%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-10 14:57:42
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-10 15:30:13
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   

function solution = read_sol(filePath)
% Open the .sol file
fileID = fopen(filePath, 'r');

% Read the file
data = textscan(fileID, '%s %f');

% Close the file
fclose(fileID);

% Extract variables and their values
varNames = data{1};
varValues = data{2};

% Store the data in a struct for easy access
solution = struct();
for i = 1:length(varNames)
    solution.(varNames{i}) = varValues(i);
end


end
