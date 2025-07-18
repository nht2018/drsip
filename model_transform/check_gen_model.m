%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2024-01-29 20:33:23
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-02-01 20:16:28
%  
%  Copyright (c) 2024, Hantao Nie, Peking University. 
%   
%   gen_model has the form of
%   min <c, x> + c0
%   st  lc <= A * x <= uc
%     F * x + g in D
%     lx <= x <= ux

% gen_model has fields 
% A, F: double matrix
% g, lc, uc, lx, ux, c: double vector
% c0: double scalar
% D: Cone

function flag = check_gen_model(model)
    % check fields
    if ~ all(isfield(model, {'A', 'F', 'g', 'lc', 'uc', 'lx', 'ux', 'c', 'c0', 'D'}))
        flag = false;
        return;
    end

    % check dimensions
    % check number of variables
    if size(model.F, 2) ~= size(model.A, 2) || length(model.lx) ~= size(model.A, 2) || length(model.ux) ~= size(model.A, 2) || length(model.c) ~= size(model.A, 2) 
        flag = false;
        return;
    end

    % check numebr of constraints 
    if size(model.A, 1) ~= length(model.lc) || size(model.A, 1) ~= length(model.uc) || size(model.F, 1) ~= length(model.g)
        flag = false;
        return;
    end

    % check cone: to be done
    flag = true;
end
    