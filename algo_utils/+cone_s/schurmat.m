%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-22 16:02:27
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function [lhs] = schurmat(lhs, H, p, model, algo)
    %% system_opt
    % 1: ldl on sparse matrix
    % 3: chol on dense matrix 
    % 4: iterative method
    At = model.At;
    K = model.K;
    
    system_opt = algo.newton.system_opt;
    if ismember(system_opt, [1, 2]) % direct method
        error('direct method not implemented yet');
        % sigma = algo.newton.sigma;
        % temp = H.coeff{p} ;
        % lhs.mat11 = schurmat_sblk(Cone.toblk(K), At, algo.par, lhs.mat11, p, temp, temp);
        % lhs.mat11 = lhs.mat11 + sigma / ( 1+ sigma) * speye(size(lhs.mat11));

        %% Warning: to use the SDPT3 function schurmat_sblk, the following code should be add to obtain algo.par
%     %% calling SDPT3 fuctions. This is preparing for using SDPT3 functions to construct schur . Only needed for SDP.
%     %%-----------------------------------------
%     %% find the combined list of non-zero
%     %% elements of Aj, j = 1:k, for each k.
%     %% IMPORTANT NOTE: Ak, C are permuted.
%     %%-----------------------------------------
%     %%
%     algo.par = struct;
%     algo.par.spdensity   = 0.4;
%     algo.par.smallblkdim = 50;
%
%     algo.par.numcolAt = length(model.b);
%     [model.At, model.c, X, Z, algo.par.permA,algo.par.permZ] = sortA(Cone.toblk(model.K), model.At, model.c, model.b, X, Z);
%     [algo.par.isspA, algo.par.nzlistA, algo.par.nzlistAsum, algo.par.isspAy, algo.par.nzlistAy] = nzlist(Cone.toblk(model.K), model.At, algo.par);
%


    elseif system_opt == 4 % iterative method
        if ~ isfield(H{p}, 'lmut')
            H{p}.lmut = @(r_) H{p}.eigs * (H{p}.coeff .* (H{p}.eigs' * r_ * H{p}.eigs)) * H{p}.eigs' + H{p}.shift .* r_;
        end
        lhs.mat11{p} = @(r_) At{p}' * mysvec(K{p}, H{p}.lmut(mysmat(K{p}, At{p} * r_)));        
    end


end