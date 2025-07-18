function addpath_sedumi
    if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64') % my mac laptop
        addpath("/Users/niehantao/Desktop/software/sedumi");
    elseif strcmp(computer, 'GLNXA64') % my linux server
        addpath("/bicmr/home/nieht/software/sedumi");
    end
end
