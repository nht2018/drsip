%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-22 15:08:12
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-09 20:38:22
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function model = sedumi2gurobi(sedumimodel)
    %% translate sedumi model to gurobi model
    if ~ isfield(sedumimodel, 'A') && isfield(sedumimodel, 'At')
        sedumimodel.A = sedumimodel.At';
    end
    assert(isfield(sedumimodel, 'K'), 'sedumi model must have field K')
    assert(isfield(sedumimodel, 'A'), 'sedumi model must have field A')
    assert(isfield(sedumimodel, 'b'), 'sedumi model must have field b')
    assert(isfield(sedumimodel, 'c'), 'sedumi model must have field c')

    model = struct();
    model.modelsense = 'min';
    model.obj = full(sedumimodel.c);
    model.A = sedumimodel.A;
    model.rhs = full(sedumimodel.b);
    model.sense = '=';
    model.lb = [zeros(sedumimodel.K.l, 1); -inf(sum(sedumimodel.K.q), 1)];
    model.ub = inf(size(sedumimodel.c));
    
    K = sedumimodel.K;
    startIndex = K.l + 1;
    n = length(sedumimodel.c);
    for i = 1:length(K.q)
        model.quadcon(i) = struct();
    end
    for i = 1:length(K.q)
        coneSize = K.q(i);
        model.quadcon(i).name = char("QC" + i);
        model.quadcon(i).q = zeros(n, 1);
        model.quadcon(i).rhs = 0;

        % The matrix part representing the quadratic terms
        idx = [startIndex:startIndex+coneSize-1]';
        % model.quadcon(i).Qc = sparse(idx, idx, [-1; ones(coneSize-1, 1)], n, n, coneSize);
        model.quadcon(i).Qrow = idx;
        model.quadcon(i).Qcol = idx;
        model.quadcon(i).Qval = [-1; ones(coneSize-1, 1)];

        model.lb(startIndex) = 0; 
        
        startIndex = startIndex + coneSize;
    end
end