%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-12 11:06:26
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-09-18 12:36:00
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   

%% compute A * H * A' (or its schur complement type transformation )for quadratic cone
% H has strcutre as 
%      H = H.diag + H.lr * diag(H.coeff) * H.lr'
%% usage 
%      lhs = linsysconstruct_cone_q(lhs, cone_size, At, H, system_opt)
%      lhs = linsysconstruct_cone_q(lhs, cone_size, At, H, system_opt, At_den, Adencol, At_sp, Aspcol)
% input 
%      lhs: struct, left hand side of the linear system
%      cone_size: cone size,
%      At: matrix
%      H: struct, must have fiels H.diag, H.lr, H.coeff
%         size(H.diag) = [sum(cone_size), 1]
%         size(H.lr) = [sum(cone_size), k], k is the rank of H.lr. Usually k = length(cone_size) or 2 * length(cone_size)
%         size(H.coeff) = [k, 1]
%         for system_opt = 'augmented'
%               we need Hinv = H.invdiag + H.invlr * diag(H.invcoeff) * H.invlr'
%               hence H must have field H.invdiag, H.invlr, H.invcoeff
%
%      system_opt: char
%            sparse: ldl on sparse matrix
%            sparse_psd: ldl on dense matrix with mat11 PSD
%            dense: chol on dense matrix 
%            augmented: augmented system
%            for 'sparse' and 'sparse_psd', At_den, Adencol, At_sp, Aspcol are needed
%      At_den: matrix, dense part of At
%      At_sp: matrix, sparse part of At
%      Adencol: int vector, dense column index of At
%      Aspcol: int vector, sparse column index of At

function [lhs] = linsysconstruct_cone_q(lhs, cone_size, At, H, system_opt, At_den, Adencol, At_sp, Aspcol)

if strcmp(system_opt, 'sparse') || strcmp(system_opt, 'sparse_psd')
    if nargin < 9
        At_den = [];
        Adencol = [];
        At_sp = At;
        Aspcol = [1: size(At, 1)];
    end
end

%% dimension check for H
% to be implemented


if strcmp(system_opt, 'sparse') 
    % low rank part
    AHlr = At' * H.lr;
    idxden = detect_dense_col(AHlr);
    idxden = intersect(idxden, find(H.coeff ~= 0));
    if ~ isempty(idxden) % separate sparse and dense handling
        idxsp = setdiff([1: size(AHlr,2)], idxden);
        AHlr_sp = AHlr(:, idxsp);
        lhs.mat11 = lhs.mat11 + AHlr_sp * spdiag(H.coeff(idxsp)) * AHlr_sp';
        lhs.mat31 = [lhs.mat31; AHlr(:, idxden)'];
        lhs.mat33_diag = [lhs.mat33_diag; -1 ./ H.coeff(idxden)];
    else % all columns are sparse
        lhs.mat11 = lhs.mat11 + AHlr * spdiag(H.coeff) * AHlr';
    end
 
    % diagonal part
    idxden = Adencol;
    idxsp = Aspcol;
    lhs.mat11 = lhs.mat11 + At_sp' * spdiag( H.diag(idxsp) ) * At_sp;
    lhs.mat21_den = [lhs.mat21_den; At_den];
    lhs.mat22_diag = [lhs.mat22_diag; - 1 ./ H.diag(idxden)];

elseif strcmp(system_opt, 'sparse_psd')
    error('not implemented yet')
elseif  strcmp(system_opt, 'dense')
    % low rank part
    AHlr = At' * H.lr;
    lhs.mat11 = lhs.mat11 + AHlr * spdiag(H.coeff) * AHlr';

    % diagonal part
    lhs.mat11 = lhs.mat11 + At' * spdiag(H.diag) * At;

elseif strcmp(system_opt, 'augmented')
    % low rank part
    lhs.mat31 = [lhs.mat31; sparse(length(cone_size), size(lhs.mat31, 2))];
    lhs.mat32 = [lhs.mat32, sparse(size(lhs.mat32, 1), size(H.invlr, 1));
                    sparse(size(H.invlr, 2), size(lhs.mat32, 2)), H.invlr'];
    lhs.mat33_diag = [lhs.mat33_diag; 1 ./ H.invcoeff];

    % diagonal part
    lhs.mat21 = [lhs.mat21; At];
    lhs.mat22_diag = [lhs.mat22_diag; - H.invdiag];
end

end