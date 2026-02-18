function plotRawRegionBranchingRatio(R, kmax, opt)
%wrapper for plotBranchingRatio
arguments
    R
    kmax
    opt.states = {'all'}
    opt.regions = [-1]
end
    profile = getRawRegionProfile(R, states=opt.states, regions=opt.regions);
    plotBranchingRatio(profile, kmax);
end

