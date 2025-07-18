%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-16 15:45:44
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-16 17:21:08
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

% probnames = ["2013_firL2L1alph"];

for probname = probnames

    if dataset == "CBLIB"
        filename = fullfile(dir_data, dataset, probname + ".mat");
        load(filename);

        fprintf("probname = %s\n", probname);
        fprintf("size(A) = [%d, %d]\n", size(model.a));
        fprintf("norm(blx(blx ~= -inf)) = %f\n", norm(model.blx(model.blx ~= -inf)));
        fprintf("norm(bux(bux ~= inf)) = %f\n", norm(model.bux(model.bux ~= inf)));
        accs_reshape = reshape(model.accs, 2, [])';
        cone_type = accs_reshape(:, 1);
        fprintf("cone type: %d\n", unique(cone_type));

        row_type = repelem(cone_type, accs_reshape(:, 2));
        row_idx1 = find(row_type == 1);
        row_idx2 = find(row_type == 2);
        row_idx3 = find(row_type == 3);
        row_idx4 = find(row_type == 4);
        row_idx5 = find(row_type == 5);
        cone_idx1 = find(cone_type == 1);
        cone_idx2 = find(cone_type == 2);
        cone_idx3 = find(cone_type == 3);
        cone_idx4 = find(cone_type == 4);
        cone_idx5 = find(cone_type == 5);


        if size(model.f, 1) >= size(model.f, 2)
            fprintf("norm(f(q part) - Identity) = %f\n", norm(model.f(row_idx4, :) - speye(length(row_idx4)), 'fro'));
        else
            fprintf("size(f, 1) < size(f, 2)\n");
        end
        % figure ;
        % spy(model.f);
        % title(probname);
        % saveas(gcf, fullfile(dir_data, dataset, probname + ".png"));
        % close gcf;

    end   


end
