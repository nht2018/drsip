%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-12 12:06:29
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-29 21:32:43
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function model = mosek2gen(mosekmodel)
    %% transform a mosek model to my general model

    % input problem is
    % min c' * x + cfix
    % s.t.  blc <= a * x <= buc
    %     f * x + g in accs
    %     blx <= x <= bux


    % tranform accs to blk
    assert(mod(length(mosekmodel.accs), 2) == 0, "length(mosekmodel.accs) must be even.");
    accs_reshape = reshape(mosekmodel.accs, 2, [])';
    cone_type = accs_reshape(:, 1);
    cone_size = accs_reshape(:, 2);
    for i = 1:length(cone_type)
        assert(cone_type(i) >= 0 && cone_type(i) <= 5, " unsupported cone type." + cone_type(i));
    end
    % cone_type = 0: free cone
    % cone_type = 1: zero cone
    % cone_type = 2: nonnegative cone
    % cone_type = 3: nonpositive cone
    % cone_type = 4: quadratic cone
    % cone_type = 5: rotated quadratic cone     
    
    row_type = repelem(cone_type, cone_size);
    row_idx0 = find(row_type == 0);
    row_idx1 = find(row_type == 1);
    row_idx2 = find(row_type == 2);
    row_idx3 = find(row_type == 3);
    row_idx4 = find(row_type == 4);
    row_idx5 = find(row_type == 5);
    cone_idx0 = find(cone_type == 0);
    cone_idx1 = find(cone_type == 1);
    cone_idx2 = find(cone_type == 2);
    cone_idx3 = find(cone_type == 3);
    cone_idx4 = find(cone_type == 4);
    cone_idx5 = find(cone_type == 5);

    % create a new my model
    model = struct();

    % move linear blocks (including zero cone) to A matrix
    model.A = [mosekmodel.a; mosekmodel.f(row_idx1, :); mosekmodel.f(row_idx2, :); mosekmodel.f(row_idx3, :)];
    model.lc = [mosekmodel.blc; - mosekmodel.g(row_idx1); - mosekmodel.g(row_idx2); -inf(length(row_idx3), 1)];
    model.uc = [mosekmodel.buc; - mosekmodel.g(row_idx1); inf(length(row_idx2), 1); - mosekmodel.g(row_idx3)];

    % concatenate the second order cones 
    model.F = mosekmodel.f([row_idx4; row_idx5], :);
    model.g = mosekmodel.g([row_idx4; row_idx5]);
    model.D = {'q', reshape(accs_reshape(cone_idx4, 2), 1, []);
                'r', reshape(accs_reshape(cone_idx5, 2), 1, [])};

    % remove empty blocks
    model.D = model.D(~cellfun(@(cone_size) sum(cone_size) == 0, model.D(:, 2)), :);
    model.D = Cone.fromblk(model.D);

    % copy the rest of the data
    model.c0 = mosekmodel.cfix;
    model.c = mosekmodel.c;
    model.lx = mosekmodel.blx;
    model.ux = mosekmodel.bux;


    
   
end