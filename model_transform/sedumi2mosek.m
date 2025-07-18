%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-22 19:10:52
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-23 12:00:09
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function model = sedumi2mosek(sedumimodel)
% sedumi2mosek - convert a sedumi model to a mosek model
    if ~ isfield(sedumimodel, 'A') && isfield(sedumimodel, 'At')
        sedumimodel.A = sedumimodel.At';
    end
    assert(isfield(sedumimodel, 'K'), 'sedumi model must have field K')
    assert(isfield(sedumimodel, 'A'), 'sedumi model must have field A')
    assert(isfield(sedumimodel, 'b'), 'sedumi model must have field b')
    assert(isfield(sedumimodel, 'c'), 'sedumi model must have field c')


    K = sedumimodel.K;
    if ~isfield(K, 'f'); K.f = 0; end
    if ~isfield(K, 'l'); K.l = 0; end
    if ~isfield(K, 'q'); K.q = []; end
    if ~isfield(K, 'r'); K.r = []; end
    if ~isfield(K, 's'); K.s = []; end

    [r, res] = mosekopt('symbcon');
    n = K.f + K.l + sum(K.q) + sum(K.r);  % number of variables exclding SDP variables
    model.c   = sedumimodel.c(1: n);
    model.a   = sedumimodel.A(:, 1: n);
    model.blc = sedumimodel.b;
    model.buc = sedumimodel.b;

    model.blx = [];
    model.bux = [];
    model.accs = [];
    model.f = sparse(0, n);
    model.g = zeros(0, 1);

    % f
    model.blx = [model.blx; - inf(K.f, 1)];
    model.bux = [model.bux; inf(K.f, 1)];

    % l
    model.blx = [model.blx; zeros(sum(K.l), 1)];
    model.bux = [model.bux; inf(sum(K.l), 1)];

    % q
    model.blx = [model.blx; -inf(sum(K.q), 1)];
    model.bux = [model.bux; inf(sum(K.q), 1)];
    temp = [repmat(res.symbcon.MSK_DOMAIN_QUADRATIC_CONE, numel(K.q), 1), K.q'];
    temp = reshape(temp', 1, []);
    model.accs = [model.accs, temp];
    model.f = [model.f; 
            sparse(sum(K.q), K.l + K.f), speye(sum(K.q)), sparse(sum(K.q), sum(K.r));];
    model.g = [model.g; zeros(sum(K.q), 1)];

    % r
    model.blx = [model.blx; -inf(sum(K.r), 1)];
    model.bux = [model.bux; inf(sum(K.r), 1)];
    temp = [repmat(res.symbcon.MSK_DOMAIN_RQUADRATIC_CONE, numel(K.r), 1), K.r'];
    temp = reshape(temp', 1, []);
    model.accs = [model.accs, temp];
    model.f = [model.f; 
            sparse(sum(K.r), K.l + K.f + sum(K.q)), speye(sum(K.r))];
    model.g = [model.g; zeros(sum(K.r), 1)];

    % s
    if ~ isempty(K.s)
        model.bardim = K.s;

        % barc
        % barc
        max_size = sum(K.s.^2);  % Maximum possible size of arrays
        model.barc.subj = zeros(max_size, 1);
        model.barc.subk = zeros(max_size, 1);
        model.barc.subl = zeros(max_size, 1);
        model.barc.val = zeros(max_size, 1);
        IndexStart = n + 1;
        index = 1;  % Index for adding elements to arrays
        for j = 1: numel(K.s)
            dim = K.s(j);
            C = reshape(sedumimodel.c(IndexStart: IndexStart + dim^2-1) , dim, dim);
            [r, c, v] = find(tril(C));
            num_elements = numel(v);
            model.barc.subj(index: index + num_elements - 1) = j;
            model.barc.subk(index: index + num_elements - 1) = r;
            model.barc.subl(index: index + num_elements - 1) = c;
            model.barc.val(index: index + num_elements - 1) = v;
            index = index + num_elements;
            IndexStart = IndexStart + dim^2;
        end
        % Trim arrays to actual size
        model.barc.subj = model.barc.subj(1: index - 1);
        model.barc.subk = model.barc.subk(1: index - 1);
        model.barc.subl = model.barc.subl(1: index - 1);
        model.barc.val = model.barc.val(1: index - 1);


        % % bara
        % model.bara.subi = zeros(max_size, 1);
        % model.bara.subj = zeros(max_size, 1);
        % model.bara.subk = zeros(max_size, 1);
        % model.bara.subl = zeros(max_size, 1);
        % model.bara.val = zeros(max_size, 1);

        
        % IndexStart = n + 1;
        % index = 1;  % Index for adding elements to arrays
        % for j = 1: numel(K.s)
        %     dim = K.s(j);
        %     for i = 1: size(sedumimodel.A, 1)
        %         A = reshape(sedumimodel.A(i, IndexStart: IndexStart + dim^2-1) , dim, dim);
        %         [r, c, v] = find(tril(A));
        %         num_elements = numel(r);
        %         model.bara.subi(index: index + num_elements - 1) = i;
        %         model.bara.subj(index: index + num_elements - 1) = j;
        %         model.bara.subk(index: index + num_elements - 1) = r;
        %         model.bara.subl(index: index + num_elements - 1) = c;
        %         model.bara.val(index: index + num_elements - 1) = v;
        %         index = index + num_elements;
        %     end
        %     IndexStart = IndexStart + dim^2;
        % end
        % % Trim arrays to actual size
        % model.bara.subi = model.bara.subi(1: index - 1);
        % model.bara.subj = model.bara.subj(1: index - 1);
        % model.bara.subk = model.bara.subk(1: index - 1);
        % model.bara.subl = model.bara.subl(1: index - 1);
        % model.bara.val = model.bara.val(1: index - 1);
    end
end





