% Set dir_data to the path where the data is stored in before running this script. There
% should be a DIMACS subdirectory under dir_data
% Set instances to be the names of the problems that you want to solve or "all"
% if you want to solve all problems in the dataset. Default is "all".


mode = "debug";
profile_on = 0;
method = "sn";
% dataset = "CBLIB";
% dataset = "DIMACS";
dataset = "NETLIB_MAT_MDO_PRSLVD";
% dataset = "MITTELMANN_MAT_MDO_PRSLVD";

dir_data = get_data_dir(dataset);
if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64')
    dir_results = "../results";
    dir_logs = "../logs";
elseif strcmp(computer, 'GLNXA64')
    dir_results = "../results_linux";
    dir_logs = "../logs_linux";
end



instances = "25FV47";

addpath_mosek;
if numel(instances) == 1 && instances == "all"
    if dataset == "DIMACS"
        probnames = ["nb", "nb_L1", "nb_L2", "nb_L2_bessel", "nql180", "nql30", "nql60", "qssp180", "qssp30", "qssp60", "sched_100_100_orig", "sched_100_100_scaled", "sched_100_50_orig", "sched_100_50_scaled",  "sched_200_100_orig", "sched_200_100_scaled", "sched_50_50_orig", "sched_50_50_scaled"];
%         probnames = ["nb", "nb_L1", "nb_L2", "nb_L2_bessel", "nql30", "nql60", "qssp30", "qssp60", "sched_100_100_orig", "sched_100_100_scaled", "sched_100_50_orig", "sched_100_50_scaled",  "sched_200_100_orig", "sched_200_100_scaled", "sched_50_50_orig", "sched_50_50_scaled"];
%         probnames = [ "sched_100_100_orig", "sched_100_100_scaled", "sched_100_50_orig", "sched_100_50_scaled", "sched_200_100_orig", "sched_200_100_scaled", "sched_50_50_orig", "sched_50_50_scaled"];
    elseif dataset == "CBLIB"
        probnames = ["2013_dsNRL", "2013_firL1", "2013_firL1Linfalph", "2013_firL1Linfeps", "2013_firL2L1alph", ...
            "2013_firL2L1eps", "2013_firL2Linfalph", "2013_firL2Linfeps", "2013_firL2a", "2013_firLinf", "2013i_wbNRL", ...
            "beam30", "beam7", "chainsing-50000-1", "chainsing-50000-2","chainsing-50000-3", "db-joint-soerensen", ...
            "db-plate-yield-line"];
       probnames = ["2013_dsNRL", "2013_firL1", "2013_firL1Linfalph", "2013_firL2L1eps","2013_firL2Linfalph", "2013i_wbNRL"] ;
       % 2013_firL2Linfeps: AHlr has dense col
       
    elseif dataset == "cbf_pre"
        probnames = ["2013_dsNRL", "2013_firL1", "2013_firL1Linfalph", "2013_firL1Linfeps", "2013_firL2L1alph", ...
            "2013_firL2L1eps", "2013_firL2Linfalph", "2013_firL2Linfeps", "2013_firL2a", "2013_firLinf", "2013i_wbNRL", ...
            "beam30", "beam7", "chainsing-50000-1", "chainsing-50000-2","chainsing-50000-3", "db-joint-soerensen", ...
            "db-plate-yield-line"];
        for i = 1: length(probnames)
            probnames(i) = probnames(i) + "pre" ;
        end
    elseif dataset == "sdplib"
        probnames = ["arch0", "arch2", "arch4", "arch8", "control1", "control10", "control11", "control2", "control3", "control4", "control5", "control6", "control7", "control8", "control9", "equalG11", "equalG51", "gpp100", "gpp124-1", "gpp124-2", "gpp124-3", "gpp124-4", "gpp250-1", "gpp250-2", "gpp250-3", "gpp250-4", "gpp500-1", "gpp500-2", "gpp500-3", "gpp500-4", "hinf1", "hinf10", "hinf11", "hinf12", "hinf13", "hinf14", "hinf15", "hinf2", "hinf3", "hinf4", "hinf5", "hinf6", "hinf7", "hinf8", "hinf9", "infd1", "infd2", "infp1", "infp2", "maxG11", "maxG32", "maxG51", "maxG55", "maxG60", "mcp100", "mcp124-1", "mcp124-2", "mcp124-3", "mcp124-4", "mcp250-1", "mcp250-2", "mcp250-3", "mcp250-4", "mcp500-1", "mcp500-2", "mcp500-3", "mcp500-4", "qap10", "qap5", "qap6", "qap7", "qap8", "qap9", "qpG11", "qpG51", "ss30", "theta1", "theta2", "theta3", "theta4", "theta5", "theta6", "thetaG11", "thetaG51", "truss1", "truss2", "truss3", "truss4", "truss5", "truss6", "truss7", "truss8"];
    elseif dataset == "NETLIB_MAT_MDO_PRSLVD"
        % remove FIT2D FIT2P GROW15 GROW22 GROW7
        probnames = strsplit("25FV47 80BAU3B ADLITTLE AFIRO AGG AGG2 AGG3 BANDM BEACONFD BLEND BNL1 BNL2 BOEING1 BOEING2 BORE3D BRANDY CAPRI CYCLE CZPROB D2Q06C D6CUBE DEGEN2 DEGEN3 DFL001 E226 ETAMACRO FFFFF800 FINNIS FIT1D FIT1P FIT2D FIT2P FORPLAN GANGES GFRD-PNC GREENBEA GREENBEB GROW15 GROW22 GROW7 ISRAEL KB2 LOTFI MAROS MAROS-R7 MODSZK1 NESM PEROLD PILOT PILOT-JA PILOT-WE PILOT4 PILOT87 PILOTNOV QAP12 QAP15 QAP8 RECIPE SC105 SC205 SC50A SC50B SCAGR25 SCAGR7 SCFXM1 SCFXM2 SCFXM3 SCORPION SCRS8 SCSD1 SCSD6 SCSD8 SCTAP1 SCTAP2 SCTAP3 SEBA SHARE1B SHARE2B SHELL SHIP04L SHIP04S SHIP08L SHIP08S SHIP12L SHIP12S SIERRA STAIR STANDATA STANDGUB STANDMPS STOCFOR1 STOCFOR2 STOCFOR3 TRUSS TUFF VTP-BASE WOOD1P WOODW");
%         probnames = strsplit(" ISRAEL KB2 LOTFI MAROS MAROS-R7 MODSZK1 NESM PEROLD PILOT PILOT-JA PILOT-WE PILOT4 PILOT87 PILOTNOV QAP12 QAP15 QAP8 RECIPE SC105 SC205 SC50A SC50B SCAGR25 SCAGR7 SCFXM1 SCFXM2 SCFXM3 SCORPION SCRS8 SCSD1 SCSD6 SCSD8 SCTAP1 SCTAP2 SCTAP3 SEBA SHARE1B SHARE2B SHELL SHIP04L SHIP04S SHIP08L SHIP08S SHIP12L SHIP12S SIERRA STAIR STANDATA STANDGUB STANDMPS STOCFOR1 STOCFOR2 STOCFOR3 TRUSS TUFF VTP-BASE WOOD1P WOODW");
    elseif dataset == "MITTELAMNN_MAT_MDO_PRSLVD"
        probnames = ["QAP15", "ex10", "neos2", "pds-40", "shs1023brazil3", "fome13", "neos3", "physiciansched3-3", "square41", "buildingenergy", "irish-electricity", "ns1644855", "rail02", "stat96v1", "chromaticindex1024-7", "l1sixm250obs", "ns1687037", "rail4284", "stat96v4", "cont1", "linf520c", "ns1688926", "s100", "stormg21000", "cont11", "neos", "nug08-3rd", "s250r10", "stp3d", "dbic1", "neos-5052403-cygnet", "nug15", "savsched1", "supportcase10", "ds-big", "neos1", "pds-100", "self-consistent", "watson2"];
    end
else
    probnames = instances;
end

%% load param file
param_file = ['params_', char(method), '.csv'];
if exist(param_file, 'file')
    Tdata = readtable(param_file, 'Delimiter', ':');
else
    Tdata = [];
end



for probname = probnames
    fprintf("*************************************************************************************\n");
    fprintf("Solving problem %s from %s dataset ...\n", probname, dataset);
    

    %% load parameters
    params = struct;
    params.tol = 1e-6;
    params.print_log = 2;
    params.max_time = 7200;
    params.iterative_refinement_max_iter = 0;
    params.print_log = 2;
    params.max_iter = 1000;
    params.mu_init = -1;
    params.AAt_solver = 'default';
    params.newton_solver = 'default';
    params.system_opt = 0;
    params.newton_tol = 1e-8;
    params.boot_method = "fom";
    params.newton_reg = 1;
    params.predcorr = 1;
%     params.rescale_option = 0;
%     params.scaling_method = 'none';
    %     params.newton_reg = 1;
    %     params.newton_solver = 'matlab_backslash';
    % params.newton_solver = 'pardiso';
%     params.init_strategy = 1;
    params.method = char(method);
    params.presolve = 'none';
    params.solve_dual = 0;
    %     params.mat_path = '/Users/niehantao/Documents/Seafile/seafile/ssnsdp_cpp/mat';
    if mode == "test"
        params.warning_on = 0;
    elseif mode == "debug"
        params.warning_on = 0;
    end
    params.log_path = '';
    % create log file if not exist and clear if exist
    mkdir(fullfile(dir_logs, dataset, probname));
    fid = fopen(fullfile(dir_logs, dataset, probname, method + ".log"), 'w');
    fclose(fid);
    diary(fullfile(dir_logs, dataset, probname, method + ".log"));
    % params.log_path = fullfile(dir_logs, dataset, "1e-" + num2str(-log(params.tol) / log(10)), probname, "drspf_" + string(params.newton_solver) + ".log");
    result_path = fullfile(dir_results, dataset, probname, method + ".txt");
    % result_path = fullfile(dir_results, dataset, "1e-" + num2str(-log(params.tol) / log(10)), probname, "drspf_" + string(params.newton_solver) + ".txt");
    % set specific params
    params = set_spec_params(Tdata, probname, params);



    if dataset == "DIMACS"
        model_orig = load(fullfile(dir_data, probname + ".mat"));
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
    elseif dataset == "CBLIB"
        if strcmp(params.presolve, 'gurobi')
            data_path = get_data_dir(dataset + "_gurobipresolve");
            if isfile(fullfile(data_path, "gen", probname + ".mat"))
                fprintf("reading %s ...\n", fullfile(data_path, "gen", probname + ".mat"));
                load(fullfile(data_path, "gen", probname + ".mat"));
            else
                error("File %s does not exist.", probname);
            end
        elseif strcmp(params.presolve, 'none')
            data_path = fullfile(dir_data);
            if isfile(fullfile(data_path, "gen", probname + ".mat"))
                fprintf("reading %s ...\n", fullfile(data_path, "gen", probname + ".mat"));
                load(fullfile(data_path, "gen", probname + ".mat"));
            elseif isfile(fullfile(data_path, "mps", probname + ".mps")) 
                fprintf("reading %s ...\n", fullfile(data_path, "mps", probname + ".mps"));
                model = gurobi_read(char(fullfile(data_path, "mps", probname + ".mps")));
                model = gurobi2gen(model);
                % filename = fullfile(data_path, "mps", probname + ".mps");
                % [r, res] = mosekopt(convertStringsToChars(strcat('read(', filename, ')')));
                % prob = res.prob;
                % model = mosek2gen(prob);
            elseif isfile(fullfile(data_path, "cbf", probname + ".cbf.gz"))
                fprintf("reading %s ...\n", fullfile(data_path, "cbf", probname + ".cbf.gz"));
                filename = fullfile(data_path, "cbf", probname + ".cbf.gz");
                [r, res] = mosekopt(convertStringsToChars(strcat('read(', filename, ')')));
                prob = res.prob;
                model = mosek2gen(prob);
            else
                error("File %s does not exist.", probname);
            end
        end
    elseif dataset == "cbf_pre"
        %         if isfile(fullfile(dir_data, dataset, "std", probname + ".mat"))
        %             load(fullfile(dir_data, dataset, "std", probname + ".mat"));
        %         else
        if isfile(fullfile(dir_data, dataset, probname + ".cbf"))
            filename = fullfile(dir_data, dataset, probname + ".cbf");
            [r, res] = mosekopt(convertStringsToChars(strcat('read(', filename, ')')));
            prob = res.prob;
            model = mosek2std(prob);
        end
    elseif dataset == "sdplib"
        if isfile(fullfile(dir_data, dataset, probname + ".dat-s"))
            [blk, At, c, b] = read_sdpa(fullfile(dir_data, dataset, probname + ".dat-s"));
            
            model = struct();
            model.b = b;
            model.K = Cone.fromblk(blk);
            model.At = MatCell(At);
            model.c = MatCell(c);
            
        else
            error("No such file: %s", fullfile(dir_data, dataset, probname + ".dat-s"));
        end
    elseif dataset == "qplib"
        model = gurobi_read(fullfile(dir_data, probname + ".lp"));
        model = gurobi2
    elseif dataset == "NETLIB_MAT_MDO_PRSLVD" || dataset == "MITTELMANN_MAT_MDO_PRSLVD"
        [A, b, c, l, u] = read_MAT_MDO_PRSLVD(fullfile(dir_data, dataset, probname));
        
        model.A = A;
        model.uc = b;
        model.lc = b;
        model.c = c;
        model.c0 = 0;
        model.g = zeros(0, 1);
        model.F = sparse(0, size(A, 2));
        model.D = {};
        model.lx = l;
        model.ux = u;
        
        model = gen2std(model) ;
        
        % idx_free = find(l == -Inf & u == Inf);
        % idx_lb = find(l ~= -Inf & u == Inf);
        % idx_ub = find(l == -Inf & u ~= Inf);
        % idx_oneside = union(idx_lb, idx_ub);
        % idx_box = find(l ~= -Inf & u ~= Inf);
        
        % model = struct();
        % model.K = cell(0, 1);
        % model.At = cell(0, 1);
        % model.c = cell(0, 1);
        % if ~isempty(idx_free)
        %     model.K{end + 1} = BasicCone('u', length(idx_free));
        %     model.At{end + 1} = A(:, idx_free)';
        %     model.c{end + 1} = c(idx_free);
        % end
        % if ~isempty(idx_box)
        %     model.K{end + 1} = BasicCone('b', length(idx_box), [l(idx_box), u(idx_box)]);
        %     model.At{end + 1} = A(:, idx_box)';
        %     model.c{end + 1} = c(idx_box);
        % end
        % if ~isempty(idx_oneside)
        %     model.K{end + 1} = BasicCone('l', length(idx_oneside));
        %     A(:, idx_ub) = -A(:, idx_ub);
        %     l(idx_ub) = -u(idx_ub);
        %     u(idx_ub) = Inf;
        %     c(idx_ub) = -c(idx_ub);
        %     b = b - A(:, idx_oneside) * l(idx_oneside);
        %     model.At{end + 1} = A(:, idx_oneside)';
        %     model.c{end + 1} = c(idx_oneside);
        % end
        % model.b = b;
        % model.c = reshape(model.c, [], 1);
        % model.At = reshape(model.At, [], 1);
    end
    
    % remove rebuntant fields in model other than At, c, b, K
    model.name = probname;
    
    
    t_drspf = tic;

    if profile_on
        profile on;
    end
    if method == "drspf"
        if mode == "test"
            try
                out  = drspf(model, params);
            catch
                continue;
            end
        elseif mode == "debug"
            [out, model, hist]  = drspf(model, params);
        end
    else 
        if mode == "test"
            try
                out  = SN(model, params);
            catch
                continue;
            end
        elseif mode == "debug"
            [out, model, hist]  = SN(model, params);
        end
    end
    diary off;
    t_drspf = toc(t_drspf);
    if profile_on
        profile off;
        profile viewer;
    end
    
    [result_dir, ~, ~] = fileparts(result_path);
    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end
    fileID = fopen(result_path, 'w');
    if strcmp(out.status, 'OPTIMAL')
        status = 's';
    else
        status = 'f';
    end
    fprintf(fileID, '%.4e\t%i\t%.4e\t%.4e\t%.4e\t%s', t_drspf, out.iter, out.gap, out.pinf, out.dinf, status);
    fclose(fileID);
end
% fprintf('%i problem(s) solved by DRSPF.\n', numel(probnames));
% maketable(dir_results, dataset, params.tol);
