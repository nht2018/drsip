function addpath_gurobi()
    if strcmp(computer, "MACI64") || strcmp(computer, "MACA64")
        addpath(genpath('/Library/gurobi1003/macos_universal2/'))
    elseif strcmp(computer, 'GLNXA64') % my linux server
        addpath(genpath('/bicmr/home/nieht/software/gurobi950/linux64'))
    end
end