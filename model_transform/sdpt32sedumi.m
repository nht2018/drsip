%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-11-02 12:05:43
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2023-11-02 12:06:04
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function model = sdpt32sedumi(sdpt3model)
    % convert a sdpt3 model to sedumi model

    K = struct;
    K.f = [];   
    K.l = [];
    K.q = [];
    K.r = [];
    K.s = [];

    A = struct;
    A.f = [];
    A.l = [];
    A.q = [];
    A.r = [];
    A.s = [];

    c = struct;
    c.f = [];
    c.l = [];
    c.q = [];
    c.r = [];
    c.s = [];



    for p = 1: size(sdpt3model.blk, 1)
        cone_type = sdpt3model.blk{p, 1};
        cone_size = sdpt3model.blk{p, 2};
        if strcmp(cone_type, 'u')
            K.f = [K.f, cone_size];
            A.f = [A.f, sdpt3model.At{p}'];
            c.f = [c.f; sdpt3model.C{p}];
        elseif strcmp(cone_type, 'l')
            K.l = [K.l, cone_size];
            A.l = [A.l, sdpt3model.At{p}'];
            c.l = [c.l; sdpt3model.C{p}];
        elseif strcmp(cone_type, 'q')
            K.q = [K.q, cone_size];
            A.q = [A.q, sdpt3model.At{p}'];
            c.q = [c.q; sdpt3model.C{p}];
        elseif strcmp(cone_type, 'r')
            K.r = [K.r, cone_size];
            A.r = [A.r, sdpt3model.At{p}'];
            c.r = [c.r; sdpt3model.C{p}];
        elseif strcmp(cone_type, 's')
            K.s = [K.s, cone_size];
            A.s = [A.s, sdpt3model.At{p}'];
            c.s = [c.s; sdpt3model.C{p}];
        end
    end

    % remove empty cone in K
    if isempty(K.f)
        K = rmfield(K, 'f');
    end
    if isempty(K.l)
        K = rmfield(K, 'l');
    end
    if isempty(K.q)
        K = rmfield(K, 'q');
    end
    if isempty(K.r)
        K = rmfield(K, 'r');
    end
    if isempty(K.s)
        K = rmfield(K, 's');
    end

    % create sedumi model
    model = struct;
    model.K = K;

    % concatenate A, c
    model.A = [A.f, A.l, A.q, A.r, A.s];
    model.c = [c.f; c.l; c.q; c.r; c.s];

    % copy b
    model.b = sdpt3model.b;

end