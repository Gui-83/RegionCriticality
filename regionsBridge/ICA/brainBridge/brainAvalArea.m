function area = brainAvalArea(brain)
    size = brainAvalSize(brain);
    lifetime = brainAvalLifetime(brain);
    area = getArea(lifetime, size);
end
