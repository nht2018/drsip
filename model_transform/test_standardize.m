dir_data = get_data_dir();

% dataset = "CBLIB_gurobipresolve";
dataset = "CBLIB";

probnames = ["2013_dsNRL", "2013_firL1", "2013_firL1Linfalph", "2013_firL1Linfeps", "2013_firL2L1alph", "2013_firL2L1eps", "2013_firL2Linfalph", "2013_firL2Linfeps", "2013_firL2a", "2013_firLinf", "2013i_wbNRL", "beam7", "beam30", "chainsing-50000-1", "chainsing-50000-2","chainsing-50000-3", "db-joint-soerensen", "db-plate-yield-line"];
% probnames = ["beam7", "beam30", "chainsing-50000-1", "chainsing-50000-2","chainsing-50000-3", "db-joint-soerensen", "db-plate-yield-line"];
probnames = ["db-joint-soerensen"];

% % clear file
% fid = fopen(dataset + ".txt", 'w');
% fclose(fid);

% diary(dataset + ".txt")

for probname = probnames
    fprintf("*************************************************************************************\n");
    fprintf("problem %s from %s dataset ...\n", probname, dataset);

    data_path = fullfile(dir_data, dataset);
    if isfile(fullfile(data_path, "gen", probname + ".mat"))
        load(fullfile(data_path, "gen", probname + ".mat"));
    elseif isfile(fullfile(data_path, "mps", probname + ".mps")) 
        model = gurobi_read(char(fullfile(data_path, "mps", probname + ".mps")));
        model = gurobi2gen(model);
        % filename = fullfile(data_path, "mps", probname + ".mps");
        % [r, res] = mosekopt(convertStringsToChars(strcat('read(', filename, ')')));
        % prob = res.prob;
        % model = mosek2gen(prob);
    elseif isfile(fullfile(data_path, "cbf", probname + ".cbf.gz"))
        filename = fullfile(data_path, "cbf", probname + ".cbf.gz");
        [r, res] = mosekopt(convertStringsToChars(strcat('read(', filename, ')')));
        prob = res.prob;
        model = mosek2gen(prob);
    else
        error("File %s does not exist.", probname);
    end
    % print_gen_model(model);
    std_model = standardize_forward(model);
    assert(check_std_model(std_model));

    fprintf("*************************************************************************************\n");
    print_std_model(std_model);

    % % clear file
    % data_path = fullfile(dir_data, dataset, "gen");

    % fid = fopen(fullfile(data_path, "info", probname + ".txt"), 'w');
    % fclose(fid);
    % diary(fullfile(data_path, "info", probname + ".txt") )
    % print_gen_model(model);
    % diary off

    % figure
    % spy(model.F)
    % saveas(gcf, fullfile(data_path, "spy", probname + "_F.png"));
    % close(gcf)

    % figure
    % spy(model.A)
    % saveas(gcf, fullfile(data_path, "spy", probname + "_A.png"));
    % close(gcf)



end
  
% diary off