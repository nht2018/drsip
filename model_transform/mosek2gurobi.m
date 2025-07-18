%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-16 14:53:30
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-22 14:57:16
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function model = mosek2gurobi(mosekmodel)
    %% transform a mosek model to gurobi model

    % input modellem is
    % min c' * x + cfix
    % s.t.  blc <= a * x <= buc
    %     f * x + g in accs
    %     blx <= x <= bux

    % output modellem is
    % min c' * x + alpha
    % s.t.  A * x = b
    %     x' * Qc * x + q' * x <= beta
    %     l <= x <= u



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

    % create a new my model
    model = struct();

    % Add slack variable y = F * x + g
    nx = size(mosekmodel.a, 2);
    ny = size(mosekmodel.f, 1);
    n = nx + ny;
    A = [mosekmodel.a, - sparse(size(mosekmodel.a, 1), ny)];
    model.obj = [mosekmodel.c; zeros(ny, 1)];
    model.A = [A; A; [mosekmodel.f, - speye(ny)]];
    model.rhs = [mosekmodel.blc; mosekmodel.blc; -mosekmodel.g];
    model.sense = [repmat('>', length(mosekmodel.buc), 1); repmat('<', length(mosekmodel.blc), 1); repmat('=', ny, 1)];

    model.lb = [mosekmodel.blx; -inf(ny, 1)];
    model.ub = [mosekmodel.bux; inf(ny, 1)];

    % conic constraints
    startIndex = nx + 1;
    n_soc = 0; % number of second order cone
    for i = 1: length(cone_type)
        coneType = cone_type(i);
        coneSize = cone_size(i);
        if coneType == 0 % free cone
            model.lb(startIndex:startIndex+coneSize-1) = -inf;
            model.ub(startIndex:startIndex+coneSize-1) = inf;
        elseif coneType == 1 % zero cone
            model.lb(startIndex:startIndex+coneSize-1) = 0;
            model.ub(startIndex:startIndex+coneSize-1) = 0;
        elseif coneType == 2 % nonnegative cone
            model.lb(startIndex:startIndex+coneSize-1) = 0;
            model.ub(startIndex:startIndex+coneSize-1) = inf;
        elseif coneType == 3 % nonpositive cone
            model.lb(startIndex:startIndex+coneSize-1) = -inf;
            model.ub(startIndex:startIndex+coneSize-1) = 0;
        elseif coneType == 4 % quadratic cone
            n_soc = n_soc + 1;
            model.quadcon(n_soc).name = char("q" + n_soc);
            model.quadcon(n_soc).Qc = sparse(n, n);
            model.quadcon(n_soc).q = zeros(n, 1);
            model.quadcon(n_soc).rhs = 0;

            % The matrn_socx part representn_socng the quadratn_socc terms
            model.quadcon(n_soc).Qc(startIndex+1:startIndex+coneSize-1, startIndex+1:startIndex+coneSize-1) = speye(coneSize-1);
            model.quadcon(n_soc).Qc(startIndex, startIndex) = -1;  

            model.lb(startIndex) = 0; 
        elseif coneType == 5 % rotated quadratic cone
            n_soc = n_soc + 1;
            model.quadcon(n_soc).name = char("r" + n_soc);
            model.quadcon(n_soc).Qc = sparse(n, n);
            model.quadcon(n_soc).q = zeros(n, 1);
            model.quadcon(n_soc).rhs = 0;

            % The matrix part representing the rotational quadratic cone
            model.quadcon(n_soc).Qc(startIndex+2:startIndex+coneSize-1, startIndex+2:startIndex+coneSize-1) = speye(coneSize-2);
            model.quadcon(n_soc).Qc(startIndex, startIndex + 1) = -1;  
            model.quadcon(n_soc).Qc(startIndex + 1, startIndex) = -1;

            model.lb(startIndex) = 0; 
            model.lb(startIndex + 1) = 0;
        else
            error("Unsupported cone type: " + coneType);
        end
        
        startIndex = startIndex + coneSize;
    end    
   
end