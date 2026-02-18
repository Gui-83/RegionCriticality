function lifetime = brainAvalLifetime(brain)
    lifetime = brain.aval_indeces(:,2) - brain.aval_indeces(:,1) + 1;
end

