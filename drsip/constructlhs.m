%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-01-09 21:03:46
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-06 22:19:10
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
function [lhs, model, algo] = constructlhs(H_rd, model, algo, params)
K = model.K;
At = model.At;
m = size(At{1}, 2);

%% initialize lhs
lhs = struct();

if algo.newton.system_opt == 1 % direct method and split dense columns
    lhs.mat11 = sparse(m, m);
    lhs.mat21 = sparse(0, m);
    lhs.mat31 = sparse(0, m);
    lhs.mat32 = []; % lhs.mat32 is always empty
    lhs.mat22_diag = []; % handle free variables
    lhs.mat33 = sparse(0, 0); % handle dense column
elseif algo.newton.system_opt == 2 % direct method and split dense rows
    assert(~ algo.model.exist_u) % to be improved for unbounded case
    dim_den = numel(algo.model.row_den) ;
    dim_sp = numel(algo.model.row_sp) ;
    lhs.mat11 = sparse(dim_den, dim_den);  % den * den
    lhs.mat12 = sparse(dim_den, dim_sp);  % den * sp
    lhs.mat2t = sparse(0, dim_sp);
    lhs.mat22_coeff = [];
elseif algo.newton.system_opt == 3 % direct method and not split dense columns 
    lhs.mat11 = sparse(m, m);
    lhs.mat21 = sparse(0, m);
    lhs.mat22_diag = []; % handle free variables
elseif algo.newton.system_opt == 4 % iterative method
    lhs.mat11 = MatCell(length(model.K)); % for iterative method, mat11 is a cell storing each block
end



%% construct the reduced system
for p = 1:length(K)
    cone = K{p};
    if strcmp(cone.type,'l') || strcmp(cone.type, 'b2l') || strcmp(cone.type, 'u2l')
        [lhs] = cone_l.schurmat(lhs, H_rd, p, model, algo);
    elseif strcmp(cone.type, 's')
        [lhs] = cone_s.schurmat(lhs, H_rd, p, model, algo);
    elseif strcmp(cone.type,'q') || strcmp(cone.type,'r2q')
        [lhs] = cone_q.schurmat(lhs, H_rd, p, model, algo);
    elseif strcmp(cone.type,'u')
        % if algo.iter == 1
        if ismember(algo.newton.system_opt, [1, 3])
            lhs.mat21 = [lhs.mat21; At{p}];
            lhs.mat22_diag = [lhs.mat22_diag; zeros(size(At{p}, 1), 1)];
        elseif algo.newton.system_opt == 4
            % lhs.mat11{p} = @(r) zeros(m);
            lhs.mat11{p} = @(r_) 0;
        end
    else
        error('Unknown cone type');
    end
end


if algo.newton.system_opt == 1
    
    lhs.mat22 = spdiag(lhs.mat22_diag);
    
    lhs.dim1 = size(lhs.mat11, 1);
    lhs.dim2 = size(lhs.mat22, 1);
    lhs.dim3 = size(lhs.mat33, 1);
    
    if strcmp(params.dense_column_strategy, 'default')
        if true
            algo.newton.dense_column_strategy = 'augmented';
        else
            algo.newton.dense_column_strategy = 'SMW';
        end
    else
        algo.newton.dense_column_strategy = params.dense_column_strategy;
    end

    if strcmp(algo.newton.dense_column_strategy, 'augmented')
        lhs.mat = [lhs.mat11, lhs.mat21', lhs.mat31';
            lhs.mat21, lhs.mat22, sparse(lhs.dim2, lhs.dim3);
            lhs.mat31, sparse(lhs.dim3, lhs.dim2), lhs.mat33]; % we need only lower part for factorization.
        
        lhs.lmut = @(x) lhs.mat * x ;
        lhs.regu_diag = [ones(lhs.dim1, 1); - ones(lhs.dim2, 1); - ones(lhs.dim3, 1)];

    else %SMW
        lhs.mat = [lhs.mat11, lhs.mat21';
        lhs.mat21, lhs.mat22]; % not including dense column
        lhs.lmut = @(x) [lhs.mat11 * x(1: lhs.dim1) - lhs.mat31' * (lhs.mat33 \ (lhs.mat31 * x(1: lhs.dim1))) + lhs.mat21' * x(lhs.dim1 + 1: lhs.dim1 + lhs.dim2);
            lhs.mat21 * x(1: lhs.dim1) + lhs.mat22 * x(lhs.dim1 + 1: lhs.dim1 + lhs.dim2)];
        lhs.regu_diag = [ones(lhs.dim1, 1); - ones(lhs.dim2, 1)];
        
    end
    
    % lhs.regu_diag = [ones(lhs.dim1, 1); sign(lhs.mat22_diag); sign(lhs.mat33_diag)];
    % lhs.regu_diag = [ones(lhs.dim1, 1); zeros(lhs.dim2, 1); zeros(lhs.dim3, 1)];
    lhs.mat_regu = spdiag(lhs.regu_diag);
    algo.newton.dim = size(lhs.mat, 1);
    %         if algo.newton.system_opt == 1 && lhs.dim2 == 0 && lhs.dim3 == 0 && strcmp(algo.newton.solver, 'ldl')
    %             fprintf("system_opt = 1, but there exists no dense column nor free variable, we use ldlchol instead of ldl\n");
    %             algo.newton.solver = 'ldlchol';
    %         end

elseif algo.newton.system_opt == 2    
    lhs.dim1 = size(lhs.mat11, 1);
    lhs.dim2 = size(lhs.mat2t, 2);
    
    lhs.lmut = @(x) [lhs.mat11 * x(1: lhs.dim1) + lhs.mat12 * x(lhs.dim1 + 1: lhs.dim1 + lhs.dim2);
    lhs.mat12' * x(1: lhs.dim1) + lhs.mat2t' * (lhs.mat22_coeff .* ( lhs.mat2t * x(lhs.dim1 + 1: lhs.dim1 + lhs.dim2) ))];
    
    [lhs.inv_mat22, lhs.mat, ~] = schur_complement(lhs.mat2t, lhs.mat22_coeff, lhs.mat11, lhs.mat12, struct(), 'augmented') ;

    lhs.regu_diag = [ones(lhs.dim1, 1)] ;
    lhs.mat_regu = spdiag(lhs.regu_diag);
    algo.newton.dim = size(lhs.mat, 1);
elseif algo.newton.system_opt == 3
    
    lhs.mat22 = spdiag(lhs.mat22_diag);
    
    lhs.dim1 = size(lhs.mat11, 1);
    lhs.dim2 = size(lhs.mat22, 1);
    lhs.dim3 = 0;
    
    lhs.mat = [lhs.mat11, lhs.mat21';
        lhs.mat21, lhs.mat22]; % we need only lower part for factorization.
    
    lhs.lmut = @(x) lhs.mat * x ;
    
    lhs.regu_diag = [ones(lhs.dim1, 1); - ones(lhs.dim2, 1)];
    lhs.mat_regu = spdiag(lhs.regu_diag);
    algo.newton.dim = size(lhs.mat, 1);
elseif algo.newton.system_opt == 4
    % idx_nu = setdiff(1:length(model.K), model.idx_u); % not bounded blocks
    % At_nu = vert_concat(model.At, idx_nu);
    if ~ isfield(algo.model, 'At_u')
        algo.model.At_u = MatCell.vert_concat(model.At, algo.model.idx_u);
    end
    lhs.dim1 = size(model.At{1}, 2);
    lhs.dim2 = sum(algo.model.size_u);
    lhs.dim3 = 0;
    lhs.mat11_lmut = @(x) MatCell.reduce_sum(apply_fn(lhs.mat11, x));
    if lhs.dim2 == 0
        lhs.mat12_lmut = @(x) zeros(lhs.dim1, 1);
        lhs.mat21_lmut = @(x) [];
    else
        lhs.mat12_lmut = @(x) algo.model.At_u' * x;
        lhs.mat21_lmut = @(x) algo.model.At_u * x;
    end
    lhs.mat22_lmut = @(x) zeros(sum(algo.model.size_u), 1);
    if lhs.dim2 == 0
        lhs.lmut = lhs.mat11_lmut;
    else
        lhs.lmut = @(x) [lhs.mat11_lmut(x(1: lhs.dim1)) + lhs.mat12_lmut(x(lhs.dim1 + 1: lhs.dim1 + lhs.dim2));
            lhs.mat21_lmut(x(1: lhs.dim1)) + lhs.mat22_lmut(x(lhs.dim1 + 1: lhs.dim1 + lhs.dim2))];
    end
    % lhs.diag = [lhs.mat11_diag; lhs.mat22_diag];
    % lhs.precond = @(x) x ./ lhs.diag;
    lhs.precond = @(x) x;
    algo.newton.dim =  lhs.dim1 + lhs.dim2;
end

end