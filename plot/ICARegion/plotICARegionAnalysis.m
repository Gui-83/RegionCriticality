function plotICARegionAnalysis(regions, opt)
arguments
  regions
  opt.significanceLevel = 0.2
  opt.dicoStep = 10
  opt.min_decade = 1
  opt.base = 1.5
  opt.threshold = 4
  opt.states = {'all'}
  opt.regions = [-1]
end
    [T, alpha, Txmin, Txmax, Tp, S, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam] = ICARegionAnalysis(regions, significanceLevel = opt.significanceLevel, dicoStep = opt.dicoStep, min_decade = opt.min_decade, base = opt.base, threshold=opt.threshold, states=opt.states, regions=opt.regions);
    plotAnalysis(T, alpha, Txmin, Txmax, Tp, S, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam);
end