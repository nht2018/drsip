function addpath_cvx()
    if strcmp(computer, "MACI64") || strcmp(computer, "MACA64")
        % rmpath(genpath('/Users/niehantao/Desktop/software/cvx/'))
        addpath(genpath('/Users/niehantao/Desktop/software/cvx3/'))
    elseif strcmp(computer, 'GLNXA64') % my linux server
        addpath(genpath('/bicmr/home/nieht/software/cvx/'))
    end
end