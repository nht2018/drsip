  
function addpath_ecos
    if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64') % my mac laptop
        addpath(genpath("/Users/niehantao/Desktop/software/ecos-matlab"))
    elseif strcmp(computer, 'GLNXA64') % my linux server
        addpath(genpath("/bicmr/home/nieht/software/ecos-matlab"))
    end
end

