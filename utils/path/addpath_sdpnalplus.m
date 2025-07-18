function addpath_sdpnalplus
    if strcmp(computer, "MACI64") || strcmp(computer, "MACA64")
        addpath(genpath('/Users/niehantao/Desktop/software/SDPNAL+v1.0'));
    elseif strcmp(computer, 'GLNXA64') % my linux server
        addpath(genpath('/bicmr/home/nieht/software/SDPNAL+v1.0'))
    end
end