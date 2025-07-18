%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-09-08 14:40:39
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-03-05 22:09:47
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   

function [d, info] = linsyssolve(lhs, rhs, info)
    %% solve linear system lhs * d = rhs
    
    assert(isfield(info, 'solver_used'));
    if strcmp(info.solver_used, 'ldlchol') 
        [d] = solve_ldlchol(rhs, info);
    elseif strcmp(info.solver_used, 'chol') 
        [d] = solve_chol(rhs, info);
    elseif strcmp(info.solver_used, 'spchol') 
        [d] = solve_spchol(rhs, info);
    elseif strcmp(info.solver_used, 'lchol') 
        [d] = solve_lchol(rhs, info) ;
    elseif strcmp(info.solver_used, 'pardiso')
        [d] = solve_pardiso(lhs, rhs, info) ;
    elseif strcmp(info.solver_used, 'ldlsparse') 
        [d] = solve_ldlsparse(rhs, info) ;
    elseif strcmp(info.solver_used, 'spldl') 
        [d] = solve_spldl(rhs, info) ;
    elseif strcmp(info.solver_used, 'ldl') 
        [d] = solve_ldl(rhs, info) ;
    elseif strcmp(info.solver_used, 'lu')
        [d] = solve_lu(rhs, info) ;
    elseif strcmp(info.solver_used, 'splu')
        [d] = solve_splu(rhs, info) ;    
    elseif strcmp(info.solver_used, 'pcg')
        assert(isfield(lhs, 'lmut'));
        assert(isfield(lhs, 'precond'));
        [d, flag, relres, iter, resvec] = pcg(lhs.lmut, rhs, info.newton_tol, 1000, lhs.precond);
        if flag ~= 0
            warning('pcg failed');
        end
        if isfield(info, 'pcg_iter')
            info.pcg_iter = info.pcg_iter + iter;
        else
            info.pcg_iter = iter;
        end
    elseif strcmp(info.solver_used, 'qmr')
        assert(isfield(lhs, 'lmut'));
        assert(isfield(lhs, 'precond'));
        [d, flag, relres, iter, resvec] = qmr(lhs.lmut, rhs, info.newton_tol, 1000, lhs.precond);
        if flag ~= 0
            warning('qmr failed');
        end
        if isfield(info, 'qmr_iter')
            info.qmr_iter = info.qmr_iter + iter;
        else
            info.qmr_iter = iter;
        end
    elseif strcmp(info.solver_used, 'matlab_backslash')
        assert(isfield(lhs, 'mat'));
        d = lhs.mat \ rhs;
    end

    % fprintf('residual: %e\n', norm(lhs.mat * d - rhs));
end