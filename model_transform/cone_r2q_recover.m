function [X, y, S] = cone_r2q_recover(X, y, S, model)
    for p = 1: length(model.K)
        cone = model.K{p};
        if strcmp(cone.type, 'r2q')
            X{p} = mexrotate_cone_r(X{p}, cone.size);
            S{p} = mexrotate_cone_r(S{p}, cone.size);
        end
    end
end

