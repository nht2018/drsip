%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-29 19:31:10
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-01-31 19:56:30
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function model = gurobi2gen(gurobimodel)
% gurobimodel is a struct with fields: A, obj, sense, rhs, lb, ub, vtype, modelname, quadcon, varnames, constrnames
    model.A = gurobimodel.A;
    model.c = gurobimodel.obj;
    model.c0 = 0;
    model.lc = zeros(length(gurobimodel.sense), 1);
    model.uc = zeros(length(gurobimodel.sense), 1);
    for i =1: length(gurobimodel.sense)
        if gurobimodel.sense(i) == '>'
            model.lc(i) = gurobimodel.rhs(i);
            model.uc(i) = inf;
        elseif gurobimodel.sense(i) == '<'
            model.uc(i) = gurobimodel.rhs(i);
            model.lc(i) = -inf;
        else
            model.lc(i) = gurobimodel.rhs(i);
            model.uc(i) = gurobimodel.rhs(i);
        end
    end

    model.lx = gurobimodel.lb;
    model.ux = gurobimodel.ub;

    cone_q_size = []; cone_q_idx = [];
    cone_r_size = []; cone_r_idx = [];


    model.F = sparse(length(gurobimodel.obj), length(gurobimodel.obj));
    for i=1: length(gurobimodel.quadcon)
        qcon = gurobimodel.quadcon(i);
        assert(qcon.sense == '<');
        assert(norm(qcon.q) < 1e-12);
        n = length(qcon.Qrow);
        assert(n == length(qcon.Qcol) && n == length(qcon.Qval));
        Q = sparse(qcon.Qrow, qcon.Qcol, qcon.Qval);
        Q = Q(qcon.Qrow, qcon.Qcol) ;
        if norm(Q - spdiag([-1; ones(n-1,1)]), 'fro' ) < 1e-12 && norm(qcon.q) < 1e-16% quadratic cone
            cone_q_idx = [cone_q_idx, i];
            cone_q_size = [cone_q_size, n];
            if abs(model.lx(qcon.Qrow(1)) - 0) < 1e-12 % romove the constrain x(1) >= 0
                model.lx(qcon.Qrow(1)) = -inf;
            end
            model.F(sub2ind(size(model.F), qcon.Qrow, qcon.Qcol)) = 1;
        elseif norm(Q - mexrotate_matrix([n]), 'fro') < 1e-12 && norm(qcon.q) < 1e-16% rotate quadratic cone
            cone_r_idx = [cone_r_idx, i];
            cone_r_size = [cone_r_size, n];
            if abs(model.lx(qcon.Qrow(1)) - 0) < 1e-12 % romove the constrain x(1), x(2) >= 0
                model.lx(qcon.Qrow(1)) = -inf;
            end
            if abs(model.lx(qcon.Qrow(2)) - 0) < 1e-12
                model.lx(qcon.Qrow(2)) = -inf;
            end
            model.F(sub2ind(size(model.F), qcon.Qrow, qcon.Qcol)) = 1;
        else
            error('Unsupported quadratic cone');
        end
    end

    % remove empty row in F
    model.F = model.F(any(model.F, 2), :);
    model.g = zeros(size(model.F, 1), 1);


    model.D = Cone() ;
    if ~isempty(cone_q_idx)
        model.D{end+1} = BasicCone('q', cone_q_size);
    end
    if ~isempty(cone_r_idx)
        model.D{end+1} = BasicCone('r', cone_r_size);
    end


    

end
        

