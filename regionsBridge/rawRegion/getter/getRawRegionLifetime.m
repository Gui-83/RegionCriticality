function T = getRawRegionLifetime(regions, opt)
arguments
    regions
    opt.states = {'all'}
    opt.regions = [-1]
end
T = [];
regionArray = regions.regions_array(:);
n = length(regionArray);
for i = 1:n
    r = regionArray(i);
    if checkState(r, opt.states) && checkRegion(r, opt.regions)
        t = r.aval_indeces(:, 2) - r.aval_indeces(:,1) + 1;;
        if isnan(t)
            continue
        end
        T = [T;t];
    end
end
end