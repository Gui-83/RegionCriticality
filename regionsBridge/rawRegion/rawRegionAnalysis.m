function [T, alpha, Txmin, Txmax, Tp, S, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam, branching, crossCorr, branchingRatiolm, autocorrDecay, mean_autocorr, autocorrlm] = rawRegionAnalysis(regions, opt)
%wrapper for analysis. It returns the results of all criticality analysis
%see README.md
arguments
  regions
  opt.significanceLevel = 0.05
  opt.dicoStep = 10
  opt.min_decade = 1
  opt.base = 2
  opt.threshold = 4
  opt.T_ref = 1/3
  opt.states = {'all'}
  opt.regions = [-1]
end
    S = getRawRegionSize(regions, states=opt.states, regions=opt.regions);
    T = getRawRegionLifetime(regions, states=opt.states, regions=opt.regions);
    A = getRawRegionArea(regions, states=opt.states, regions=opt.regions);
    
    timeDependendentSize = getRawRegionSizeTimeDependent(regions, states=opt.states, regions=opt.regions);
    [ST, lengthST] = separateAvalSizeTimeDependent(timeDependendentSize);
    profile = getRawRegionProfile(regions, states=opt.states, regions=opt.regions);

    [alpha, Txmin, Txmax, Tp, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam, branching, crossCorr, branchingRatiolm, autocorrDecay, mean_autocorr, autocorrlm] = analysis(S, T, A, ST, lengthST, profile, significanceLevel=opt.significanceLevel, dicoStep=opt.dicoStep, min_decade=opt.min_decade, base=opt.base, threshold=opt.threshold, T_ref=opt.T_ref);
end