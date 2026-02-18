function A = getICARegionArea(regions, opt)
arguments
    regions
    opt.states = {'all'}
    opt.regions = [-1]
end
    T = getICARegionLifetime(regions, states=opt.states, regions=opt.regions);
    S = getICARegionSize(regions, states=opt.states, regions=opt.regions);
    A = getArea(T, S);
end