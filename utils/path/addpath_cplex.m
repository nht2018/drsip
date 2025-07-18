function addpath_cplex()
    if strcmp(computer, "MACI64") || strcmp(computer, "MACA64")
        addpath("/Applications/CPLEX_Studio_Community2211/cplex/bin/x86-64_osx/")
    elseif strcmp(computer, 'GLNXA64') % my linux server
        addpath("/bicmr/home/nieht/software/CPLEX/cplex/bin/x86-64_linux")
    end
end