%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-16 16:59:29
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-29 12:23:54
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%  

%We use mosek to read CBLIB data. 
% The mosek model is formulated as
% min c'*x
% st blc <= A * x <= buc
%    blx <= x <= bux
%    f * x + g \in accs
% We observe that data in CBLIB has similar structure as follows:
% 1. lower part of f approximately equals identity matrix
% 2. lower part of g is a vector with all zeros
% 3. A, blc, buc are all empty. Instead, accs begins with zeros cone. that is, mosek move linear constraints to the upper part of conic constraints.
% 4. all blx = -inf, bux = inf




% Set dir_data to the path where the data is storede in before running this script. There
% should be a DIMACS subdirectory under dir_data
% Set instances to be the names of the problems that you want to solve or "all"
% if you want to solve all problems in the dataset. Default is "all".


dir_data = '../../data/sdp_data';

dataset = "CBLIB";
instances = "all";   
if numel(instances) == 1 && instances == "all"
    if dataset == "CBLIB"
        % remove beam30 and beam 7
        probnames = ["2013_dsNRL", "2013_firL1", "2013_firL1Linfalph", "2013_firL1Linfeps", "2013_firL2L1alph", ...
            "2013_firL2L1eps", "2013_firL2Linfalph", "2013_firL2Linfeps", "2013_firL2a", "2013_firLinf", "2013i_wbNRL", ...
            "chainsing-50000-1", "chainsing-50000-2","chainsing-50000-3", "db-joint-soerensen", ...
            "db-plate-yield-line"];
    end
else
    probnames = instances;
end



for probname = probnames
    fprintf("*************************************************************************************\n");
    fprintf("Presolving problem %s from %s dataset ...\n", probname, dataset);

    if dataset == "CBLIB"
        if isfile(fullfile(dir_data, dataset, probname + ".mat"))
            load(fullfile(dir_data, dataset, probname + ".mat"));
            model = mosek2std(model);
        elseif isfile(fullfile(dir_data, dataset, probname + ".cbf.gz"))
            filename = fullfile(dir_data, dataset, probname + ".cbf.gz");
            [r, res] = mosekopt(convertStringsToChars(strcat('read(', filename, ')')));
            prob = res.prob;
            model = mosek2std(prob);
        else
            error("File %s does not exist.", probname);
        end
    end

%     rmfield(model, 'At_mat');
%     rmfield(model, 'c_mat');
%     rmfield(model, 'blk');

    print_std_model(model);

    % save(fullfile(dir_data, dataset, "std", probname + ".mat"), 'model', '-v7.3');
    try
        save(fullfile(dir_data, dataset, "std", probname + ".mat"), 'model');
    catch
        print("Saving problem %s failed. try to save with -v7.3 option.", probname);
        save(fullfile(dir_data, dataset, "std", probname + ".mat"), 'model', '-v7.3');
    end
        


end
  