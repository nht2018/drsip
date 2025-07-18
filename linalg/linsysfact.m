%%
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-10-31 23:14:52
%  Description:
%
%  Copyright (c) 2023, Hantao Nie, Peking University.
%

function [info, flag] = linsysfact(lhs, info)
%% factorize the linear system lhs * d = rhs
assert(ismember(info.solver, {'auto', 'pardiso', 'pcg', 'qmr', 'chol', 'spchol', 'lchol', 'ldlchol', 'ldlsparse', 'spldl', 'ldl', 'splu', 'lu', 'matlab_backslash'}), 'unkown solver %s', info.solver);
assert(isfield(info, 'solver'));
info.solver_used = info.solver;
if ~isfield(lhs, 'mat_regu')
    lhs.mat_regu = [];
end
if strcmp(info.solver_used, 'auto')
    % we consider cholesky and ldl for PSD case and non-PSD case respectively
    if all(diag(lhs.mat) > 0) % PSD case
        if nnz(lhs.mat) / numel(lhs.mat) > 0.1 % dense case
            info.solver_used = 'chol';
        else % sparse case
            if exist('ldlchol', 'file') % use SuiteSparse ldlchol
                info.solver_used = 'ldlchol';
            elseif exist('lchol', 'file') % use SuiteSparse lchol
                info.solver_used = 'lchol';
            else
                info.solver_used = 'chol'; % use matlab chol
            end
        end
    else % non-PSD case
        if nnz(lhs.mat) / numel(lhs.mat) > 0.5
            lhs.mat = full(lhs.mat);
            info.solver_used = 'ldl';
        else
            lhs.mat = sparse(lhs.mat);
            info.solver_used = 'spldl';
            % info.solver_used = 'ldlsparse'; % seems some bud in ldlsparse I don't know why
        end
    
    end
end


if strcmp(info.solver_used, 'pardiso')
    [info, flag] = fact_pardiso(lhs.mat, info, lhs.mat_regu);  
elseif strcmp(info.solver_used, 'pcg')
    info.pcg_iter = 0;
    % info.R = ichol(lhs.mat); % preconditioner
    % info.R = ichol(lhs.mat, struct('michol','on'));
elseif strcmp(info.solver_used, 'qmr')
    info.qmr_iter = 0;
    % info.R = ichol(lhs.mat); % preconditioner
    % info.R = ichol(lhs.mat, struct('michol','on'));  
elseif strcmp(info.solver_used, 'ldlchol')  % ld for positive definite matrix
    [info, flag] = fact_ldlchol(lhs.mat, info, lhs.mat_regu);
    if flag ~= 0
        warning('ldlchol fail, switch to splu');
        info.solver_used = 'splu'; 
    end
elseif strcmp(info.solver_used, 'chol')
    [info, flag] = fact_chol(lhs.mat, info, lhs.mat_regu);
    if flag ~= 0
        warning('chol fail, switch to lu');
        info.solver_used = 'lu'; 
    end
elseif strcmp(info.solver_used, 'spchol')
    [info, flag] = fact_spchol(lhs.mat, info, lhs.mat_regu);
    if flag ~= 0
        warning('spchol fail, switch to splu'); 
        info.solver_used = 'splu'; 
    end
elseif strcmp(info.solver_used, 'lchol')
    [info, flag] = fact_lchol(lhs.mat, info, lhs.mat_regu);
    if flag ~= 0
        warning('lchol fail, switch to splu');
        info.solver_used = 'splu'; 
    end
elseif strcmp(info.solver_used, 'ldlsparse')  % ld for not neccessarily positive definite matrix
    [info, flag] = fact_ldlsparse(lhs.mat, info, lhs.mat_regu);
    if flag ~= 0
        warning('ldlsparse fail, switch to splu'); 
        info.solver_used = 'splu'; 
    end
elseif strcmp(info.solver_used, 'spldl')
    [info, flag] = fact_spldl(lhs.mat, info, lhs.mat_regu);
    if flag ~= 0
        warning('ldl fail, switch to splu');
        info.solver_used = 'splu'; 
    end
elseif strcmp(info.solver_used, 'ldl')
    [info, flag] = fact_ldl(lhs.mat, info, lhs.mat_regu);
    if flag ~= 0
        warning('ldl fail, switch to lu');
        info.solver_used = 'lu'; 
    end
end


% Once Cholesky or LDL fail, we switch to LU factorization in all later iterations
if strcmp(info.solver_used, 'lu')
    [info, flag] = fact_lu(lhs.mat, info, lhs.mat_regu);
    if flag ~= 0
        warning('lu fail, switch to matlab backslash');
        info.solver_used = 'matlab_backslash'; 
    end
elseif strcmp(info.solver_used, 'splu')
    [info, flag] = fact_splu(lhs.mat, info, lhs.mat_regu);
    if flag ~= 0
        warning('splu fail, switch to matlab backslash');
        info.solver_used = 'matlab_backslash'; 
    end
end

% if LU still fail, we switch to matlab backslash 
if strcmp(info.solver_used, 'matlab_backslash')

end


end

