%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-16 21:34:19
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-12-17 23:50:21
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   

addpath_gurobi;
dataset = "qplib";
dir_data = '/Users/niehantao/Desktop/software/data/sdp_data/qplib/html/lp';
dir_out = '/Users/niehantao/Desktop/software/data/sdp_data/qplib_grbpsv';
if ~exist(dir_out, 'dir')
    mkdir(dir_out);
end
probnames = ["QPLIB_2456", "QPLIB_2468", "QPLIB_2482", "QPLIB_2519", "QPLIB_2626", "QPLIB_2676", "QPLIB_2784", "QPLIB_2862", "QPLIB_3029", "QPLIB_3088", "QPLIB_3105", "QPLIB_3185", "QPLIB_3312", "QPLIB_8495", "QPLIB_8500", "QPLIB_8515", "QPLIB_8547", "QPLIB_8559", "QPLIB_8567", "QPLIB_8602", "QPLIB_8616", "QPLIB_8785", "QPLIB_8790", "QPLIB_8792", "QPLIB_8845", "QPLIB_8906", "QPLIB_8938", "QPLIB_8991", "QPLIB_9002", "QPLIB_9008", "QPLIB_10034", "QPLIB_10038", "QPLIB_8585", "QPLIB_8803"] ;

% probnames = ["QPLIB_2456", "QPLIB_2468", "QPLIB_2482", "QPLIB_2519", "QPLIB_2626", "QPLIB_2676", "QPLIB_2784", "QPLIB_2862", "QPLIB_3029", "QPLIB_3088", "QPLIB_3105", "QPLIB_3185", "QPLIB_3312", "QPLIB_8495", "QPLIB_8500", "QPLIB_8515", "QPLIB_8547", "QPLIB_8559", "QPLIB_8567", "QPLIB_8602", "QPLIB_8616", "QPLIB_8785", "QPLIB_8790", "QPLIB_8792", "QPLIB_8845", "QPLIB_8906", "QPLIB_8938", "QPLIB_8991", "QPLIB_9002", "QPLIB_9008", "QPLIB_10034", "QPLIB_10038", "QPLIB_QP8585", "QPLIB_QP8803", "QPLIB_bdry3", "QPLIB_co5_1", "QPLIB_cont5", "QPLIB_dist3", "QPLIB_dtoc", "QPLIB_hub1", "QPLIB_rqp1", "QPLIB_twod"];

probnames = ["QPLIB_2456"] ;

for p = 1: length(probnames)
    probname = probnames(p);
    model = gurobi_read(char(fullfile(dir_data, probname + ".lp")));

    if ismember(probname, ["QPLIB_8585", "QPLIB_8803"] )
        params = struct();
        params.NonConvex = 2;
        model_psv = gurobi_presolve(model, params);
    else
        model_psv = gurobi_presolve(model);
    end
    gurobi_write(model_psv, char(fullfile(dir_out, probname + ".lp")));
end

