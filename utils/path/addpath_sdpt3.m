  function addpath_sdpt3
    if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64')  % my mac laptop
        SDPT3Home = '/Users/niehantao/Desktop/software/SDPT3-4.0';
        eval(['addpath ',strcat(SDPT3Home,'/')]);
        eval(['addpath ',strcat(SDPT3Home,'/Solver')]);
        eval(['addpath ',strcat(SDPT3Home,'/Solver/Mexfun')]);
        eval(['addpath ',strcat(SDPT3Home,'/HSDSolver')]);
        eval(['addpath ',strcat(SDPT3Home,'/Examples')]);
    elseif strcmp(computer, 'GLNXA64') % my linux server
        SDPT3Home = '/bicmr/home/nieht/software/SDPT3-4.0';
        eval(['addpath ',strcat(SDPT3Home,'/')]);
        eval(['addpath ',strcat(SDPT3Home,'/Solver')]);
        eval(['addpath ',strcat(SDPT3Home,'/Solver/Mexfun')]);
        eval(['addpath ',strcat(SDPT3Home,'/HSDSolver')]);
        eval(['addpath ',strcat(SDPT3Home,'/Examples')]);
    end
end
