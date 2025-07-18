function addpath_scs
    if strcmp(computer, "MACI64") || strcmp(computer, "MACA64")
        addpath(genpath('/Users/niehantao/Desktop/software/scs-matlab/'))
        addpath(genpath('/Users/niehantao/Desktop/software/scs/'))
    elseif strcmp(computer, 'GLNXA64') % my linux server
        addpath(genpath('/bicmr/home/nieht/software/scs-matlab/'))
        addpath(genpath('/bicmr/home/nieht/software/scs/'))
    end
end