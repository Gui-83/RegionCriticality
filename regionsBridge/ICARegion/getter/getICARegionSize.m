function S = getICARegionSize(regions, opt)
arguments
    regions
    opt.states = {'all'}
    opt.regions = [-1]
end
S = [];
regionArray = regions.regions_array(:);
n = length(regionArray);
for i = 1:n
    r = regionArray(i);
    if checkState(r, opt.states) && checkRegion(r, opt.regions)
        s = r.ICaval_sizes;
        if isnan(s)
            continue
        end
        S = [S;s];
    end
end
end