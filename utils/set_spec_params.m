function new_params = set_spec_params(Tdata, probname, params)

    new_params = params;
    
    % skip when data is empty
    if isempty(Tdata)
        return;
    end
    
    % find valid row
    par = Tdata(strcmp(Tdata.PROBLEM, probname), :);
    
    if isempty(par)
        return;
    end
    
    % set param
    val = par.PARAMS;
    val = split(val{1}, ',');
    
    if numel(val) == 1
        return;
    end

    for i = 1:2:numel(val)
        % val{i+1} has format num1;num2;num3
        val_c = split(val{i+1}, ';');
        new_params.(val{i}) = str2double(val_c);
        if isnan(new_params.(val{i}))
            new_params.(val{i}) = char(val_c);
        end
    end
    
end
    