  
function addpath_suitesparse
    if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64') % my mac laptop
        run("/Users/niehantao/Desktop/software/SuiteSparse/SuiteSparse_paths")
    elseif strcmp(computer, 'GLNXA64') % my linux server
        run("/bicmr/home/nieht/software/SuiteSparse/SuiteSparse_paths");
    end
end

