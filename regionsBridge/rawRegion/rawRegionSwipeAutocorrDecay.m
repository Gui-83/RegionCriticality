function chi = rawRegionSwipeAutocorrDecay(regions, opt)
%swipe the reference time parameter T_ref to find the best fit of slow
%autocorelation decay. See Tian et .al Theoretical foundations of studying
%criticality in the brain

%This is a wrapper for swipeAutocorrelationDecay, it just extract data and
%feed it to the function
arguments
    regions
    opt.step = 10
    opt.rSquaredMin = 0.9
    opt.interpn = 1000
    opt.states = {'all'}
    opt.regions = [-1]
end
    timeDependendentSize = getRawRegionSizeTimeDependent(regions, states=opt.states, regions=opt.regions);
    ST = separateAvalSizeTimeDependent(timeDependendentSize);
    chi = swipeAutocorrelationDecay(ST, step=opt.step, interpn=opt.interpn, rSquaredMin=opt.rSquaredMin);
end

