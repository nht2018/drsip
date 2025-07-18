mex_folder = 'mexfun';
files = dir(fullfile(mex_folder, '*.c'));
files = {files.name};

for i = 1:length(files)
    filename = files{i};
    command = ['mex -R2018a -outdir ', mex_folder, ' CXX_FLAGS="-Xclang -fopenmp" LDFLAGS="$LDFLAGS -lomp" CXXOPTIMFLAGS="$CXXOPTIMFLAGS -Xclang -fopenmp -O3" -I/usr/local/include ', fullfile(mex_folder, filename)];
    disp(command);
    eval(command);   
end