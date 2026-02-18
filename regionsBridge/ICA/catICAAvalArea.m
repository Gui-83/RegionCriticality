function area = catICAAvalArea(region)
    lifetime = catICAAvalLifetime(region);
    size = catICAAvalSize(region);
    area = getArea(lifetime, size);
end

