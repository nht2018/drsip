%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-10-29 21:26:22
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-01-29 20:46:39
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function gen_model = std2gen(std_model)
%% transform a standard model to a general model
% Input:
%   std_model: standard model
%   min <C, X> 
%   st   A * x == b
%        x in K
%              std_model has fields At, blk, C, b
% Output:
%   gen_model: general model
%   min <c, x> + c0
%   st  lc <= A * x <= uc
%     F * x + g in D
%     lx <= x <= ux
%              gen_model has fields A, F, g, D, lc, uc, lx, ux, C, c0
%              gen_model can be feed to drs_sqlp

    At = MatCell.vert_concat(std_model.At);
    [n, m] = size(At);
    gen_model = struct;
    if isfield(std_model, 'name')
        gen_model.name = std_model.name;
    end
    gen_model.A = At';
    gen_model.lc = std_model.b;
    gen_model.uc = std_model.b;


    gen_model.lx = [];
    gen_model.ux = [];
    gen_model.D = {};
    ind = 1;
    gen_model.F = sparse(0, n);
    for p = 1: size(std_model.blk)
        len = sum(std_model.blk{p, 2});
        if strcmp(std_model.blk{p, 1}, 'l')
            gen_model.lx = [gen_model.lx; zeros(len, 1)];
            gen_model.ux = [gen_model.ux; inf(len , 1)];
        elseif strcmp(std_model.blk{p, 1}, 'u')
            gen_model.lx = [gen_model.lx; -inf(len, 1)];
            gen_model.ux = [gen_model.ux; inf(len, 1)];
        elseif strcmp(std_model.blk{p, 1}, 'q')
            gen_model.lx = [gen_model.lx; -inf(len, 1)];
            gen_model.ux = [gen_model.ux; inf(len, 1)];
            gen_model.D = {gen_model.D;
                                BasicCone('q', std_model.blk{p, 2})};
            gen_model.F = [gen_model.F;
                            sparse(len, ind - 1), speye(len, len), sparse(len, n - ind - len + 1)];
        elseif strcmp(std_model.blk{p, 1}, 'r')
            gen_model.lx = [gen_model.lx; -inf(len, 1)];
            gen_model.ux = [gen_model.ux; inf(len, 1)];
            gen_model.D = {gen_model.D; 
                                BasicCone('r', std_model.blk{p, 2})};
            gen_model.F = [gen_model.F;
                            sparse(len, ind - 1), speye(len, len), sparse(len, n - ind - len + 1)];
        else
            error('support l, u, q, r blocks only');
        end
        ind = ind + len;
    end

    gen_model.g = zeros(size(gen_model.F, 1), 1);
    gen_model.c = MatCell.vert_concat(std_model.c);
    gen_model.c0 = 0;

    
end