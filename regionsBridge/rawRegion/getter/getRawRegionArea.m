function A = getRawRegionArea(regions, opt)
arguments
    regions
    opt.states = {'all'}
    opt.regions = [-1]
end
    T = getRawRegionLifetime(regions, states=opt.states, regions=opt.regions);
    S = getRawRegionSize(regions, states=opt.states, regions=opt.regions);
    A = getArea(T, S);
end

