%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-16 15:45:44
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-16 16:40:02
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
% transform cbf my model and save as .mat file


addpath_mosek;% we need mosek to read cbf files

dir_data = "../../data/sdp_data";
dataset = "CBLIB";


probnames = ["2013_dsNRL", "2013_firL1", "2013_firL1Linfalph", "2013_firL1Linfeps", "2013_firL2L1alph", ...
    "2013_firL2L1eps", "2013_firL2Linfalph", "2013_firL2Linfeps", "2013_firL2a", "2013_firLinf", "2013i_wbNRL", ...
    "beam30", "beam7", "chainsing-50000-1", "chainsing-50000-2","chainsing-50000-3", "db-joint-soerensen", ...
    "db-plate-yield-line"];


for probname = probnames

    if dataset == "CBLIB"
        filename = fullfile(dir_data, dataset, probname + ".cbf.gz");
        fprintf("Reading %s\n", filename);
        if isfile(filename)
            [r, res] = mosekopt(convertStringsToChars(strcat('read(', filename, ')')));
        else
            error("File %s not found.\n", filename);
        end
        model = res.prob;
        save(fullfile(dir_data, dataset, probname + ".mat"), 'model');

    end   


end
