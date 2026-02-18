function chi = plotRawRegionSwipeAutocorrDecay(regions, opt)
%wrapper for rawRegionSwipeAutocorrDecay, it plots the slow decay of
%autocorelation parameters vs reference time.
arguments
    regions
    opt.step = 10
    opt.rSquaredMin = 0.9
    opt.interpn = 1000
    opt.states = {'all'}
    opt.regions = [-1]
end
    chi = rawRegionSwipeAutocorrDecay(regions, step=opt.step, rSquaredMin=opt.rSquaredMin, interpn=opt.interpn, states=opt.states, regions=opt.regions);
    x = linspace(0, 1, opt.step);
    plot(x, chi);
    ylabel('chi');
    xlabel('t');
end

