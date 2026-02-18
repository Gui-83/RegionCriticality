function plotRawRegionAnalysis(regions, opt)
%wrapper for rawRegionAnalysis and plotAnalysis (see README.md)
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
    [T, alpha, Txmin, Txmax, Tp, S, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam, branching, crossCorr, branchingRatiolm, autocorrDecay, mean_autocorr, autocorrlm] = rawRegionAnalysis(regions, significanceLevel = opt.significanceLevel, dicoStep = opt.dicoStep, min_decade = opt.min_decade, base = opt.base, threshold=opt.threshold, T_ref=opt.T_ref, states=opt.states, regions=opt.regions);
    plotAnalysis(T, alpha, Txmin, Txmax, Tp, S, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam, branching, crossCorr, branchingRatiolm, autocorrDecay, mean_autocorr, autocorrlm, opt.T_ref);
end
