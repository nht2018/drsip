%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-16 21:34:19
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-17 23:53:22
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   

addpath_gurobi;
dataset = "qplib";
dir_data = '/Users/niehantao/Desktop/software/data/sdp_data/qplib/html/lp';
dir_data_psv = '/Users/niehantao/Desktop/software/data/sdp_data/qplib_grbpsv';
probnames = ["QPLIB_2456", "QPLIB_2468", "QPLIB_2482", "QPLIB_2519", "QPLIB_2626", "QPLIB_2676", "QPLIB_2784", "QPLIB_2862", "QPLIB_3029", "QPLIB_3088", "QPLIB_3105", "QPLIB_3185", "QPLIB_3312", "QPLIB_8495", "QPLIB_8500", "QPLIB_8515", "QPLIB_8547", "QPLIB_8559", "QPLIB_8567", "QPLIB_8602", "QPLIB_8616", "QPLIB_8785", "QPLIB_8790", "QPLIB_8792", "QPLIB_8845", "QPLIB_8906", "QPLIB_8938", "QPLIB_8991", "QPLIB_9002", "QPLIB_9008", "QPLIB_10034", "QPLIB_10038", "QPLIB_8585", "QPLIB_8803"];

% probnames = ["QPLIB_2456", "QPLIB_2468", "QPLIB_2482", "QPLIB_2519", "QPLIB_2626", "QPLIB_2676", "QPLIB_2784", "QPLIB_2862", "QPLIB_3029", "QPLIB_3088", "QPLIB_3105", "QPLIB_3185", "QPLIB_3312", "QPLIB_8495", "QPLIB_8500", "QPLIB_8515", "QPLIB_8547", "QPLIB_8559", "QPLIB_8567", "QPLIB_8602", "QPLIB_8616", "QPLIB_8785", "QPLIB_8790", "QPLIB_8792", "QPLIB_8845", "QPLIB_8906", "QPLIB_8938", "QPLIB_8991", "QPLIB_9002", "QPLIB_9008", "QPLIB_10034", "QPLIB_10038", "QPLIB_QP8585", "QPLIB_QP8803", "QPLIB_bdry3", "QPLIB_co5_1", "QPLIB_cont5", "QPLIB_dist3", "QPLIB_dtoc", "QPLIB_hub1", "QPLIB_rqp1", "QPLIB_twod"];

% probnames = ["QPLIB_8803"] ;

fid = fopen('qplib.txt', 'w');
for p = 1: length(probnames)
    %% read original data
    probname = probnames(p);
    model1 = gurobi_read(char(fullfile(dir_data, probname + ".lp")));
    if isfield(model1, 'Q') && isfield(model1, 'quadcon')
        type = 'qcqp';
    elseif isfield(model1, 'Q')
        type = 'qp';
    elseif isfield(model1, 'quadcon')
        type = 'qcp';
    else
        type = 'lp';
    end
    info1 = {'%12s', 'prob', probname ;
            '%5s', 'type', type;
            '%7s', 'var', num2str(size(model1.A, 2));
            '%7s', 'lb', num2str(sum(model1.lb ~= -inf & model1.ub == inf));
            '%7s', 'ub', num2str(sum(model1.lb == -inf & model1.ub ~= inf));
            '%7s', 'box', num2str(sum(model1.lb ~= -inf & model1.ub ~= inf));
            '%7s', 'free', num2str(sum(model1.lb == -inf & model1.ub == inf));
            '%7s', 'con', num2str(size(model1.A, 1));
            '%7s', 'nnz', num2str(nnz(model1.A));
    };

    %% read presolved data
    model2 = gurobi_read(char(fullfile(dir_data_psv, probname + ".lp")));
    if isfield(model2, 'Q') && isfield(model2, 'quadcon')
        type = 'qcqp';
    elseif isfield(model2, 'Q')
        type = 'qp';
    elseif isfield(model2, 'quadcon')
        type = 'qcp';
    else
        type = 'lp';
    end
    info2 = {
            % '%12s', 'prob', probname ;
            '%5s', 'type', type;
            '%7s', 'var', num2str(size(model2.A, 2));
            '%7s', 'lb', num2str(sum(model2.lb ~= -inf & model2.ub == inf));
            '%7s', 'ub', num2str(sum(model2.lb == -inf & model2.ub ~= inf));
            '%7s', 'box', num2str(sum(model2.lb ~= -inf & model2.ub ~= inf));
            '%7s', 'free', num2str(sum(model2.lb == -inf & model2.ub == inf));
            '%7s', 'con', num2str(size(model2.A, 1));
            '%7s', 'nnz', num2str(nnz(model2.A));
    };


    info3 = {
            '%7s', 'var', num2str(size(model2.A, 2) - size(model1.A, 2));
            '%7s', 'con', num2str(size(model2.A, 1) - size(model1.A, 1));
            '%7s', 'nnz', num2str(nnz(model2.A) - nnz(model1.A));
    };
    % print info
    if p == 1
        for i = 1:size(info1, 1)
            fprintf(fid, info1{i, 1}, info1{i, 2});
        end
        fprintf(fid, " | ");
        for i = 1:size(info2, 1)
            fprintf(fid, info2{i, 1}, info2{i, 2});
        end
        fprintf(fid, " | ");
        for i = 1:size(info3, 1)
            fprintf(fid, info3{i, 1}, info3{i, 2});
        end
        fprintf(fid, '\n');

    end

    for i = 1:size(info1, 1)
        fprintf(fid, info1{i, 1}, info1{i, 3});
    end
    fprintf(fid, " | ");
    for i = 1:size(info2, 1)
        fprintf(fid, info2{i, 1}, info2{i, 3});
    end
    fprintf(fid, " | ");
    for i = 1:size(info3, 1)
        fprintf(fid, info3{i, 1}, info3{i, 3});
    end
    fprintf(fid, '\n');
end

fclose(fid);