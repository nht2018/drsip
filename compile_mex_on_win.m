mex_folder = 'mexfun';
files = dir(fullfile(mex_folder, '*.c'));
files = {files.name};

for i = 1:length(files)
    filename = files{i};
    command = ['mex -R2018a -outdir ', mex_folder, ' ', fullfile(mex_folder, filename), ' COMPFLAGS="$COMPFLAGS /O3 /openmp" '];
    disp(command);
    eval(command);   
end