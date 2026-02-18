function timeDependendentSize = getICARegionSizeTimeDependent(regions, opt)
arguments
    regions
    opt.states = {'all'}
    opt.regions = [-1]
end
timeDependendentSize = [];
regionArray = regions.regions_array(:);
n = length(regionArray);
for i = 1:n
    r = regionArray(i);
    if checkState(r, opt.states) && checkRegion(r, opt.regions)
        t = r.ICaval_timeDependendentSize;
        if isnan(t)
            continue
        end
        if t(1) ~= 0
                t = [0;t];
        end
        timeDependendentSize = [timeDependendentSize; t];
    end
end
end