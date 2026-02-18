function lifetime = catICAAvalLifetime(region)
    lifetime = [];
    for i = 1:numel(region.brain_array)
        brain = region.brain_array(i);
        brain_lifetime = brainAvalLifetime(brain);
        lifetime = [lifetime; brain_lifetime];
    end
end

