%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-02-26 11:26:20
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-26 11:29:53
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-16 16:59:29
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-29 12:23:54
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%  




% Set dir_data to the path where the data is storede in before running this script. There
% should be a DIMACS subdirectory under dir_data
% Set instances to be the names of the problems that you want to solve or "all"
% if you want to solve all problems in the dataset. Default is "all".



dataset = "DIMACS";
dir_data = get_data_dir();

probnames = ["nb", "nb_L1", "nb_L2", "nb_L2_bessel", "nql180", "nql30", "nql60", "qssp180", "qssp30", "qssp60", "sched_100_100_orig", "sched_100_100_scaled", "sched_100_50_orig", "sched_100_50_scaled",  "sched_200_100_orig", "sched_200_100_scaled", "sched_50_50_orig", "sched_50_50_scaled"];
% probnames = ["nb"];

% % clear file
% fid = fopen(dataset + ".txt", 'w');
% fclose(fid);

% diary(dataset + ".txt")

for probname = probnames
    fprintf("*************************************************************************************\n");
    fprintf("problem %s from %s dataset ...\n", probname, dataset);

    data_path = fullfile(dir_data, dataset);
    model_orig = load(fullfile(data_path, probname + ".mat"));
    model_orig.c = reshape(model_orig.c, [], 1);
    model_orig.b = reshape(model_orig.b, [], 1);
    if isfield(model_orig, 'A')
        model_orig.At = model_orig.A';
    else
        model_orig.A = model_orig.At';
    end
    
    K = cell(2, 1);
    At = cell(2, 1);
    c = cell(2, 1);
    K{1} = BasicCone('l', model_orig.K.l);
    K{2} = BasicCone('q', model_orig.K.q);
    c{1} = model_orig.c(1 : model_orig.K.l);
    c{2} = model_orig.c(model_orig.K.l + 1 : end);
    At{1} = model_orig.At(1 : model_orig.K.l, :);
    At{2} = model_orig.At(model_orig.K.l + 1 : end, :);
    b = model_orig.b;
    
    model = struct;
    model.At = MatCell(At);
    model.c = MatCell(c);
    model.b = b;
    model.K = Cone(K);
    clear At K b c

%     rmfield(model, 'At_mat');
%     rmfield(model, 'c_mat');
%     rmfield(model, 'blk');

    % clear file
    data_path = fullfile(dir_data, dataset);

    fid = fopen(fullfile(data_path, "info", probname + ".txt"), 'w');
    fclose(fid);
    diary(fullfile(data_path, "info", probname + ".txt") )
    print_std_model(model);
    diary off

    % figure
    % spy(model.F)
    % saveas(gcf, fullfile(data_path, "spy", probname + "_F.png"));
    % close(gcf)

    figure
    spy(model_orig.A)
    saveas(gcf, fullfile(data_path, "spy", probname + "_A.png"));
    close(gcf)

    % data_path = fullfile(dir_data, dataset, "gen");

    % % save(fullfile(dir_data, dataset, "std", probname + ".mat"), 'model', '-v7.3');
    % try
    %     save(fullfile(data_path, probname + ".mat"), 'model');
    % catch
    %     fprintf("Saving problem %s failed. try to save with -v7.3 option.\n", probname);
    %     save(fullfile(data_path, probname + ".mat"), 'model', '-v7.3');
    % end
        


end
  
% diary off