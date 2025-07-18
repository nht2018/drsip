%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-12 12:22:02
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-22 14:58:36
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function mymodel = ecos2std(ecosmodel)
    %% transform a mosek model to standand model

    % input problem is
    % min c' * x 
    % s.t.  A * x == b
    %     h - G * x in K



    % create a new my model
    model = struct();

    % move linear blocks (including zero cone) to A matrix
    model.A = ecosmodel.A
    model.lc = [ecosmodel.blc; - ecosmodel.g(row_idx1); - ecosmodel.g(row_idx2); -inf(length(row_idx3), 1)];
    model.uc = [ecosmodel.buc; - ecosmodel.g(row_idx1); inf(length(row_idx2), 1); - ecosmodel.g(row_idx3)];

    % concatenate the second order cones 
    model.F = ecosmodel.f([row_idx4; row_idx5], :);
    model.g = ecosmodel.g([row_idx4; row_idx5]);
    model.D = {'q', reshape(accs_reshape(cone_idx4, 2), 1, []);
                'r', reshape(accs_reshape(cone_idx5, 2), 1, [])};

    % remove empty blocks
    model.D = model.D(~cellfun(@(cone_size) sum(cone_size) == 0, model.D(:, 2)), :);

    % copy the rest of the data
    model.c0 = ecosmodel.cfix;
    model.c = ecosmodel.c;
    model.lx = ecosmodel.blx;
    model.ux = ecosmodel.bux;

    % gen2std my model to be standard form
    std_model = gen2std(model);


    mymodel = struct();
    mymodel.At_mat = std_model.At_int;
    mymodel.c_mat = std_model.c_int;
    mymodel.b = std_model.b_int;
    mymodel.blk = std_model.blk_int;

    % transform At, C to cell
    mymodel.At = mat2cell(mymodel.At_mat, cellfun(@(cone_size) sum(cone_size), mymodel.blk(:, 2)), size(mymodel.At_mat, 2));
    mymodel.C = mat2cell(mymodel.c_mat, cellfun(@(cone_size) sum(cone_size), mymodel.blk(:, 2)), 1);
    mymodel.K = Cone.fromblk(mymodel.blk);
 

    
   
end