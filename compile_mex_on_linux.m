%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-19 14:22:25
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
mex_folder = 'mexfun';
files = dir(fullfile(mex_folder, '*.c'));
files = {files.name};
for i = 1:length(files)
    filename = files{i};
    filename = filename(1:end-2); % remove .c
    % command1 = ['mex', ' -R2018a -c CFLAGS="-fopenmp -fPIC -O3" ', ' -outdir ', mex_folder, ' ', fullfile(mex_folder, filename), '.c'];
    % disp(command1);
    % eval(command1);

    % command2 = ['mex', ' -R2018a', ' -L/bicmr/soft/MATLAB/R2022b/sys/os/glnxa64 -liomp5 -lpthread ', ' -outdir ', mex_folder, ' ', fullfile(mex_folder, filename), '.o'];
    % disp(command2);
    % eval(command2);

    command3 = ['mex', ' -R2018a CFLAGS="-fopenmp -fPIC -O3"  -L/bicmr/soft/MATLAB/R2022b/sys/os/glnxa64 -liomp5 -lpthread ', ' -outdir ', mex_folder, ' ', fullfile(mex_folder, filename), '.c'];
    disp(command3);
    eval(command3);
end
