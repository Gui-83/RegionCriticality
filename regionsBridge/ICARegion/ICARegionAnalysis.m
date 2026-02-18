function [T, alpha, Txmin, Txmax, Tp, S, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam] = ICARegionAnalysis(regions, opt)
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
    S = getICARegionSize(regions, states=opt.states, regions=opt.regions);
    T = getICARegionLifetime(regions, states=opt.states, regions=opt.regions);
    A = getICARegionArea(regions, states=opt.states, regions=opt.regions);
    
    timeDependendentSize = getICARegionSizeTimeDependent(regions, states=opt.states, regions=opt.regions);
    ST = separateAvalSizeTimeDependent(timeDependendentSize);
    ST_AVG = AvalAverageSizeTimeDependent(ST);

    [alpha, Txmin, Txmax, Tp, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam] = analysis(S, T, A, ST_AVG, significanceLevel=opt.significanceLevel, dicoStep=opt.dicoStep, min_decade=opt.min_decade, base=opt.base, threshold=opt.threshold);
end