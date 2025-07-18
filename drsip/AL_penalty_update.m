%%  
%  Author: Hantao Nie (nht@pku.edu.cn)
%  Date: 2023-11-08 20:06:56
%  LastEditors: Hantao Nie (nht@pku.edu.cn)
%  LastEditTime: 2024-01-04 10:15:57
%  Description: 
%  
%  Copyright (c) 2023, Hantao Nie, Peking University. 
%   
function algo = AL_penalty_update(hist, algo, params)
    smean = @geo_mean;
    iter = algo.iter;
    
    sitr = iter - params.AL_penalty_update_iter + 1;
    sitr = max(1, sitr);
    avg_pvd = smean(hist.pvd(sitr:iter));
    avg_dvp = smean(hist.dvp(sitr:iter));
    
    if avg_dvp > params.AL_penalty_update_threshold
        algo.AL_penalty = algo.AL_penalty * (params.AL_penalty_update_factor);
    elseif avg_pvd > params.AL_penalty_update_threshold
        algo.AL_penalty = algo.AL_penalty / (params.AL_penalty_update_factor);
    else
        
    end
    
    algo.AL_penalty = min(params.AL_penalty_max, max(params.AL_penalty_min, algo.AL_penalty));
end