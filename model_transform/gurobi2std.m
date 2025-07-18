function model = gurobi2std(gurobimodel)
    model = gurobi2gen(gurobimodel);
    model = gen2std(model);
end