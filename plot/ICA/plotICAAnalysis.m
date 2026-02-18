function plotICAAnalysis(region, opt)
arguments
  region
  opt.significanceLevel = 0.2
  opt.dicoStep = 10
  opt.min_decade = 1
  opt.base = 1.5
  opt.threshold = 4
end
    [T, alpha, Txmin, Txmax, Tp, S, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam] = ICAAnalysis(region, significanceLevel = opt.significanceLevel, dicoStep = opt.dicoStep, min_decade = opt.min_decade, base = opt.base, threshold=opt.threshold);
    plotAnalysis(T, alpha, Txmin, Txmax, Tp, S, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam);
end