%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-12-23 16:58:52
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-01-31 20:35:33
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-12 16:35:51
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function print_gen_model(model, fid)

    if nargin < 2
        fid = 1;
    end

    assert(all(isfield(model, {'A', 'F', 'D', 'g', 'lc', 'uc', 'lx', 'ux', 'c', 'c0'})), ...
        "gen_model must have fields 'A', 'F', 'D', 'g', 'lc', 'uc', 'lx', 'ux', 'c', 'c0'");


    fprintf(fid, "[Constrants]\n");
    A = model.A;
    fprintf(fid, "\t[Matrix A] size: %i x %i, norm = %.4e, nnz = %.4e, sparsity = %.4e\n", ...
        size(A, 1), size(A, 2), norm(A, 'fro'), nnz(A), nnz(A) / numel(A));
    fprintf(fid, "\tequations: %i\n", sum(model.lc == model.uc));
    fprintf(fid, "\tone side: %i\n", sum(model.lc ~= model.uc & ( model.lc == -inf | model.uc == inf)));
    fprintf(fid, "\ttwo sides: %i\n", sum(model.lc ~= model.uc & ( model.lc ~= -inf & model.uc ~= inf)));
    

    fprintf(fid, "[Conic Constraints]\n");
    F = model.F;
    if size(F, 1) == size(F, 2) && norm(F - speye(size(F)), 'fro') < 1e-16
        fprintf(fid, "\t[Matrix F] is identity\n");
    else
        fprintf(fid, "\t[Matrix F] size: %i x %i, norm = %.4e, nnz = %.4e, sparsity = %.4e\n", ...
            size(F, 1), size(F, 2), norm(F, 'fro'), nnz(F), nnz(F) / numel(F));
        if isdiag(F)
            fprintf(fid, "\tF is diagonal\n");
        end
    end

    if norm(model.g) == 0
        fprintf(fid, "\tg is all zero\n");
    else
        fprintf(fid, "norm(g) = %.4e\n", norm(model.g));
    end

    D = model.D;
    fprintf(fid, "[Cone D]\n");
    print_tol = 10;
    for p = 1: length(D)
        cone = D{p};
        fprintf(fid, "\tblock %i: '%s' ", p, cone.type);
        if length(cone.size) <= print_tol  % print all if not too many
            fprintf(fid, "[ ");
            for j = 1: length(cone.size)
                fprintf(fid, "%i ", cone.size(j));
            end
            fprintf(fid, "]\n");
        else    % print first and last few if too many
            fprintf(fid, "[ ");
            n_head = ceil(print_tol / 2); n_tail = print_tol - n_head;
            for j = 1: n_head
                fprintf(fid, "%i ", cone.size(j));
            end
            fprintf(fid, "... ");
            for j = length(cone.size) - n_tail + 1: length(cone.size)
                fprintf(fid, "%i ", cone.size(j));
            end
            fprintf(fid, "]\n");
            fprintf(fid, "\t\t(Cartesian product of %i cones with sizes ranging from %i to %i)\n", ...
                numel(cone.size), min(cone.size), max(cone.size));
        end  
    end 

    fprintf(fid, "[Variable Bounds]\n");
    fprintf(fid, "\tfixed: %i\n", sum(model.lx == model.ux));
    fprintf(fid, "\tfree: %i\n", sum(model.lx == -inf & model.ux == inf));
    fprintf(fid, "\tlower: %i\n", sum(model.lx ~= -inf & model.ux == inf));
    fprintf(fid, "\tupper: %i\n", sum(model.lx == -inf & model.ux ~= inf));
    fprintf(fid, "\tbox: %i\n", sum(model.lx ~= -inf & model.ux ~= inf & model.lx < model.ux));


    
end