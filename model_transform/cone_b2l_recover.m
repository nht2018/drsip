function [X, y, S, dobj_shift] = cone_b2l_recover(X, y, S, model)
    dobj_shift = 0;
    for p = 1: length(model.K)
        cone = model.K{p};
        if strcmp(cone.type, 'b2l')
            dobj_shift = dobj_shift + sum(S{p}(1: cone.size / 2) .* cone.params(:, 1)) - sum(S{p}(cone.size / 2 + 1: cone.size) .* cone.params(:, 2));
            X{p} = X{p}(1: cone.size / 2) + cone.params(:, 1);
            % S{p} = S{p}(1: cone.size / 2) - S{p}(cone.size / 2 + 1: cone.size);
        end
    end
    y = y(1: end - model.n_box);
end