function T = getICARegionLifetime(regions, opt)
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
        if isempty(r.ICaval_indeces)
            continue
        end
        t = r.ICaval_indeces(:, 2) - r.ICaval_indeces(:,1) + 1;
        if isnan(t)
            continue
        end
        T = [T;t];
    end
end
end