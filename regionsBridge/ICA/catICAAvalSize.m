function size = catICAAvalSize(region)
    size = [];
    for i = 1:numel(region.brain_array)
        brain = region.brain_array(i);
        brain_size = brainAvalSize(brain);
        size = [size; brain_size];
    end
end

