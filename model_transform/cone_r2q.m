function model = cone_r2q(model)
% cone_r2q - convert a rotated cone to a rotated quadratic cone
%
%   model = cone_r2cone_q(model)
%

K = model.K;
for p = 1: length(K)
    cone = K{p};
    if cone.type == 'r'
        rotmat = mexrotate_matrix(cone.size);
        model.At{p} = rotmat * model.At{p};
        model.c{p} = rotmat * model.c{p};
        model.K{p}.type = 'r2q';
    end
end