function [stdmodel, transform] = gen2std(model)
    %% gen2std a general socp 

    % input problem is
    % min c' * x + c0
    % s.t.  lc <= A * x <= uc
    %     F * x + g in D
    %     lx <= x <= ux
    %

    % output problem is
    % min c' * x + c0
    % s.t. A * x = b
    %       x \in K



    t0 = tic;
    transform = struct;

    assert(check_gen_model(model));
    
    assert(all(model.lc <= model.uc));
    % remove rows with lc = -inf, uc = inf (i.e., unconstrained rows)
    unconstrained_rows = find(model.lc == -inf & model.uc == inf);
    model.A(unconstrained_rows, :) = [];
    model.lc(unconstrained_rows) = [];
    model.uc(unconstrained_rows) = [];

    % remove cols with lx == ux (i.e., fixed variables)
    fixed_cols = find(model.lx == model.ux);
    model.g = model.g + model.F(:, fixed_cols) * model.lx(fixed_cols);
    model.F(:, fixed_cols) = [];
    model.c0 = model.c0 + sum(model.lx(fixed_cols) .* model.c(fixed_cols));
    model.c(fixed_cols) = [];
    model.A(:, fixed_cols) = [];
    model.lx(fixed_cols) = [];
    model.ux(fixed_cols) = [];
    
    %% classify constraints and variables
    con_idx.e = find(model.lc == model.uc);                   % equation constraints
    con_idx.ie = find(model.lc ~= model.uc);                  % inequality constraints
    con_idx.l = find( (model.lc ~= -inf & model.uc == inf) ); % lower bound constraints
    con_idx.u = find( (model.lc == -inf & model.uc ~= inf) ); % upper bound constraints
    con_idx.b = find( (model.lc ~= -inf & model.uc ~= inf & model.lc ~= model.uc) ); % boxed constraints

    n_con.total = size(model.lc, 1);
    n_con.e = length(con_idx.e);    
    n_con.ie = length(con_idx.ie); 
    n_con.l = length(con_idx.l);
    n_con.u = length(con_idx.u);
    n_con.b = length(con_idx.b);

    var_idx.l = find( (model.lx ~= -inf & model.ux == inf) ); % lower bound variables
    var_idx.u = find( (model.lx == -inf & model.ux ~= inf) ); % upper bound variables
    var_idx.b = find( (model.lx ~= -inf & model.ux ~= inf) ); % boxed variables
    var_idx.f = find( (model.lx == -inf & model.ux == inf) ); % free variables
    var_idx.nf = find( ~ (model.lx == -inf & model.ux == inf) ); % not free variables


    n_var.total = size(model.lx, 1);
    n_var.l = length(var_idx.l);
    n_var.u = length(var_idx.u);
    n_var.b = length(var_idx.b);
    n_var.f = length(var_idx.f);
    n_var.nf = length(var_idx.nf);
    n_var.soc = size(model.F, 1);

    %% Introduce slack variables
    %compute transformation matrix for constraints and variables
    transform.Lambda_w = ones(n_con.total, 1);
    transform.Lambda_w(con_idx.u) = -1;
    transform.Lambda_w_ie = transform.Lambda_w(con_idx.ie);
    transform.Lambda_x = ones(n_var.total, 1);
    transform.Lambda_x(var_idx.u) = -1;
    transform.h_w = model.lc;
    transform.h_w(con_idx.u) = model.uc(con_idx.u);
    transform.h_x = model.lx;
    transform.h_x(var_idx.u) = model.ux(var_idx.u);
    transform.h_x(var_idx.f) = 0;

    %gen2std constraints and variables
    model.c = transform.Lambda_x .* model.c;
    model.c0 = model.c0 + transform.h_x' * model.c;
    model.b = transform.Lambda_w .* (transform.h_w - model.A * transform.h_x);
    model.A = spdiag(transform.Lambda_w) * model.A * spdiag(transform.Lambda_x);
    model.g = model.g + model.F * transform.h_x;
    model.F = model.F * spdiag(transform.Lambda_x); 


    % Here we obtain the following gen2stdd form SOCP
    % min c' * x + c0
    % s.t. A * x - w = b
    %      F * x - y = g
    %      x(var_idx_box) + x_bar(var_idx_box) = ux(var_idx_box) - lx(var_idx_box)
    %      w(con_idx_box) + w_bar(con_idx_box) = uc(con_idx_box) - lc(con_idx_box)
    %      w >= 0, w_bar(con_idx_box) >= 0
    %      x(var_idx_nf) >= 0, x_bar(var_idx_box) >= 0, x(var_idx_f) free
    %      y in D


    %% Utilize the structure of F
    % find digonal subblock in F
    F_pattern = (model.F ~= 0);
    nnz_per_row = full(sum(F_pattern, 2));
    cadidate_row = find(nnz_per_row == 1);
    nnz_per_col = sum(F_pattern, 1);
    cadidate_col = intersect(find(nnz_per_col == 1), var_idx.f);

    Omega = reshape(cadidate_col, length(cadidate_col), 1);
    Phi = zeros(length(cadidate_col), 1);
    for i = 1: length(cadidate_col)
        Phi(i) = find(F_pattern(:, cadidate_col(i)));
    end
    clear candiate_row candiate_col nnz_per_row nnz_per_col F_pattern


    Phi_bar = setdiff([1: size(model.F, 1)]', Phi);
    Omega_bar = setdiff(var_idx.f, Omega);

    %% eliminate
    model.F_Phi_Omega = diag(model.F(Phi, Omega));
    model.A(:, Omega) = model.A(:, Omega) * spdiag(1./ model.F_Phi_Omega);
    model.F = model.F(Phi_bar, :);
    model.F(:, Omega) = 0;

    %% Omega is the free variables
    var_idx.f = Omega_bar;
    var_idx.q = Omega;


    %% rearrange varaibles
    model.A = [model.A(:, var_idx.nf), model.A(:, var_idx.f), model.A(:, var_idx.q)];
    model.F = [model.F(:, var_idx.nf), model.F(:, var_idx.f), model.F(:, var_idx.q)];
    model.lx = [model.lx(var_idx.nf); model.lx(var_idx.f); model.lx(var_idx.q)];
    model.ux = [model.ux(var_idx.nf); model.ux(var_idx.f); model.ux(var_idx.q)];
    model.c = [model.c(var_idx.nf); model.c(var_idx.f); model.c(var_idx.q)];

    var_idx.nf = [1: length(var_idx.nf)]';
    var_idx.f = [length(var_idx.nf) + 1: length(var_idx.nf) + length(var_idx.f)]';
    var_idx.q = [length(var_idx.nf) + length(var_idx.f) + 1: length(var_idx.nf) + length(var_idx.f) + length(var_idx.q)]';
    var_idx.l = find( (model.lx ~= -inf & model.ux == inf) ); % lower bound variables
    var_idx.u = find( (model.lx == -inf & model.ux ~= inf) ); % upper bound variables
    var_idx.b = find( (model.lx ~= -inf & model.ux ~= inf) ); % boxed variables


    n_var.q = length(var_idx.q);


    % some rows in F are eliminated, hence recompute n_var.soc
    assert(n_var.q + size(model.F, 1) == n_var.soc)
    n_var.soc = size(model.F, 1);

    assert(size(model.A, 2) == n_var.total && size(model.A, 1) == n_con.total);
    assert(size(model.F, 2) == n_var.total );



    %% concatenate A and F
    model.Abar = [model.A; 
                   model.F];


    % for csc format , use At is faster for computing A * D * At in newton system
    model.Atbar = model.Abar';

    % fprintf("[after stadandize] Abar norm: %e\n", norm(model.Abar, 'fro'));


    %% compute constraint index in inequility part
    % model.lc_ie = model.lc(con_idx.ie);
    % model.uc_ie = model.uc(con_idx.ie);
    % con_ie_idx.l = find( (model.lc_ie ~= -inf & model.uc_ie == inf) ); % lower bound constraints
    % con_ie_idx.u = find( (model.lc_ie == -inf & model.uc_ie ~= inf) ); % upper bound constraints
    % con_ie_idx.b = find( (model.lc_ie ~= -inf & model.uc_ie ~= inf & model.lc_ie ~= model.uc_ie) ); % boxed constraints

    %% build the internal representation of the problem

    temp1 = sparse(n_con.total, n_con.ie);
    temp1(con_idx.ie, :) = speye(n_con.ie); % temp1: [I_e: 0 ; I_e^c: I] , size = [n_con.total, n_con.ie]
    temp1 = [temp1; sparse(n_var.soc, n_con.ie)] ; 
    temp2 = sparse(n_var.b, n_var.total);  
    temp2(:, var_idx.b) = speye(n_var.b);  
    % temp2: [\Gamma_b: I , \Gamma_l \cup \Gamma_u : 0, \Gamma_f: 0 ], size = [n_var.b, n_var.total]
    temp3 = sparse(n_con.b, n_con.total);
    temp3(:, con_idx.b) = speye(n_con.b); 
    temp3 = temp3(:, con_idx.ie);  %temp3: [I_b: I, I_b^c: 0], size = [n_con.b, n_con.ie]

    % "int" stand for "internal" here
    % size of A_int is
    % [size(Abar, 1) + size(F, 1) + n_var.b + n_con.b, size(Abar, 2) + n_var.soc + n_con.ie + n_var.b + n_con.b]
    % size of At_int is
    % [size(Abar, 2) + n_var.soc + n_con.ie + n_var.b + n_con.b, size(Abar, 1) + size(F, 1) + n_var.b + n_con.b]
    % where size(Abar) = [n_con.total + n_var.soc, n_var.total]



    model.At_int = [model.Atbar, temp2', sparse(n_var.total, n_con.b);
                    sparse(n_var.soc, n_con.total), - speye(n_var.soc), sparse(n_var.soc, n_var.b), sparse(n_var.soc, n_con.b);
                    -temp1', sparse(n_con.ie, n_var.b), temp3';
                    sparse(n_var.b, n_con.total + n_var.soc), speye(n_var.b), sparse(n_var.b, n_con.b);
                    sparse(n_con.b, n_con.total + n_var.soc + n_var.b), speye(n_con.b)];
                    
    model.A_int = [model.Abar, [sparse(n_con.total, n_var.soc); - speye(n_var.soc)], - temp1, sparse(n_con.total + n_var.soc, n_var.b+ n_con.b);
                temp2, sparse(n_var.b, n_var.soc), sparse(n_var.b, n_con.ie), speye(n_var.b), sparse(n_var.b, n_con.b);
                sparse(n_con.b, n_var.total), sparse(n_con.b, n_var.soc), temp3, sparse(n_con.b, n_var.b), speye(n_con.b)] ;
    
    assert(norm(model.A_int' - model.At_int, 'fro') < 1e-14);
    
    model.b_int = [model.b + model.A(:, var_idx.q) * model.g(Phi, :);
                        - model.g(Phi_bar); 
                        model.ux(var_idx.b) - model.lx(var_idx.b); 
                        model.uc(con_idx.b) - model.lc(con_idx.b)];
    c_y = zeros(n_var.soc + n_var.q, 1);
    c_y(Phi) = model.c(var_idx.q) ./ model.F_Phi_Omega;
    model.c0_int = model.c0 - sum(c_y(Phi) .* model.g(Phi));
    model.c_int = [model.c(var_idx.nf);
                   model.c(var_idx.f);
                    c_y;
                    zeros(n_con.ie, 1);
                    zeros(n_var.b, 1);
                    zeros(n_con.b, 1)];

    model.blk_int = [{'l', n_var.nf;    % x_nf
                  'u', length(var_idx.f)};     % x_f
                   Cone.toblk(model.D);           % y
                  {'l', n_con.ie;   % w
                  'l', n_var.b;    % x_bar(var_idx_box)
                  'l', n_con.b}];   % w_bar(con_idx_box)

    
    % remove empty blocks
    model.blk_int = model.blk_int(~cellfun(@(cone_size) sum(cone_size) == 0, model.blk_int(:, 2)), :);

    
    % model.int_idx_f = [n_var.nf + 1 : n_var.nf + length(Phi_bar)]'; % index of free variables in internal model, i.e, Omega_bar
    % model.int_idx_q = [n_var.nf + length(Phi_bar) + 1 : n_var.nf + length(Phi_bar) + n_var.soc]'; % index of quadratic variables in internal model
    % model.int_idx_l =   [[1: n_var.nf], ...
    %             [n_var.nf + length(Phi_bar) + n_var.soc + 1 : n_var.nf + length(Phi_bar) + n_var.soc + n_con.ie]]'; % index of linear variables in internal model



        
    % transform At, C to cell
    model1 = struct;
    model1.K = Cone.fromblk(model.blk_int);
    block_size = cellfun(@(cone) sum(cone), model.blk_int(:, 2));
    model1.At = MatCell.vert_split(model.At_int, block_size);
    model1.c = MatCell.vert_split(model.c_int, block_size);
    model1.b = model.b_int;


    % concatenate the blocks of same type
    stdmodel = struct;
    [blk_type, ~, blk_idx] = unique(model.blk_int(:, 1));
    blk = cell(length(blk_type), 2);
    stdmodel.At = cell(length(blk_type), 1);
    stdmodel.c = cell(length(blk_type), 1);
    for i = 1: length(blk_type)
        blk{i, 1} = blk_type{i};
        blk{i, 2} = horzcat(model.blk_int{blk_idx == i, 2});
        if strcmp(blk_type{i}, 'l') || strcmp(blk_type{i}, 'u')
            blk{i, 2} = sum(blk{i, 2});
        end
        stdmodel.At{i} = vertcat(model1.At{blk_idx == i});
        stdmodel.c{i} = vertcat(model1.c{blk_idx == i});
    end
    stdmodel.K = Cone.fromblk(blk);
    stdmodel.b = model1.b;

    transform.time = toc(t0);
end