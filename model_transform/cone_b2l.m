
function [model] = cone_b2l(model)
n_box = 0;
model.c_int = model.c;
model.At_int = model.At;
model.b_int = model.b;
model.c0 = 0;
model.At_box = [];
for p = 1: length(model.K)
    cone = model.K{p};
    if strcmp(cone.type, 'b')
        n_box = n_box + cone.size;

        lb = cone.params(:, 1);
        ub = cone.params(:, 2);
        
        assert(all(lb > -inf) && all(ub < inf) && all(lb < ub));
        model.K{p} = BasicCone('b2l', cone.size * 2, cone.params);
        model.b = model.b - model.At{p}' * lb;
        model.c0 = model.c0 + model.c{p}' * lb;
        model.At_box = [model.At_box; model.At{p}];
        % We don't modify A here but modify it in the funtion AXfun and Atyfun
        for q = 1: length(model.K)
            if q ~= p
                % if issparse(model.At{p})
                model.At_int{q} = [model.At{q}, sparse(size(model.At{q}, 1), cone.size)];
                % else
                %    model.At{q} = [model.At{q}, zeros(size(model.At{q}, 1), cone.size)];
                % end
            else
                % if issparse(model.At{p})
                model.At_int{q} = [model.At{q}, speye(cone.size, cone.size);
                    sparse(cone.size, size(model.At{q}, 2) ), speye(cone.size, cone.size)];
                % else
                %     model.At{q} = [model.At{q}, eye(cone.size, cone.size);
                %                        zeros(cone.size, size(model.At{q}, 2)), eye(cone.size, cone.size)];
                % end
            end
        end
        
        if issparse(model.c{p})
            model.c_int{p} = [model.c{p}; sparse(size(model.c{p}))];
        else
            model.c_int{p} = [model.c{p}; zeros(size(model.c{p}))];
        end
        
        model.b_int = [model.b; model.K{p}.params(:, 2) - model.K{p}.params(:, 1)] ;
    end
end
model.n_box = n_box;
end