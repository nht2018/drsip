function addpath_superscs
    if strcmp(computer, "MACI64") || strcmp(computer, "MACA64")
        addpath(genpath('/Users/niehantao/Desktop/software/superscs/'))
    elseif strcmp(computer, 'GLNXA64') % my linux server
        addpath(genpath('/bicmr/home/nieht/software/superscs/'))
    end
end