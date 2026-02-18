function [alpha, Txmin, Txmax, Tp, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam, branching, crossCorr, branchingRatiolm, autocorrDecay, mean_autocorr, autocorrlm] = analysis(S, T, A, ST, lengthST, profile, opt)
arguments
  S %size
  T %lifetime
  A %area
  ST %time dependent size (see separateAvalSizeTimeDependent function)
  lengthST %length of each element of ST (see separateAvalSizeTimeDependent function)
  profile %profile of the spking activity (see region class)
  opt.significanceLevel = 0.95 %significance level for the Kolmogorov-Smirnov test used in the power law fit
  opt.dicoStep = 10 % number of step to perform the dichotomy in the power law fit (see docs/powerLawFit.md for more information)
  opt.min_decade = 1 % minimum number of decade between xmin and xmax for the power law fit (see docs/powerLawFit.md for more information)
  opt.base = 1.5 % base of the grid search for the power law fit (see docs/powerLawFit.md for more information)
  opt.threshold = 4 %minimun lifetime of avalanches to be considered in the collapse shape computation
  opt.maxLag = 5 %maximum lag for the branching ration computation
  opt.T_ref = 1/3; %reference time for the slow decay of autocorrelation
  opt.verbose = true
end
    if isempty(T)
        [alpha, Txmin, Txmax, Tp, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam, branching, crossCorr, branchingRatiolm, autocorrDecay, mean_autocorr, autocorrlm] = deal(NaN);
        return
    end
    [alpha, Txmin, Txmax, Tp] = swipeMLEDiscretePowerLawBounds(T, opt.significanceLevel, opt.dicoStep, opt.min_decade, opt.base, opt.verbose);
    [beta, Sxmin, Sxmax, Sp] = swipeMLEDiscretePowerLawBounds(S, opt.significanceLevel, opt.dicoStep, opt.min_decade, opt.base, opt.verbose);
    Alm = fitPowerFunction(A);
    gam = (alpha - 1)/(beta - 1);
    
    ST_AVG = AvalAverageSizeTimeDependent(ST, lengthST);
    [x_shape, shape, T_shape] = transformCollapseShape(ST_AVG);
    [T_shape, shape] = lifetimeThresholdCollapseShapeTransformed(T_shape, shape, opt.threshold);
    shape_gam = fitCollapseShape(T_shape, shape, opt.dicoStep);
  
    [branching, crossCorr, branchingRatiolm] = branchingRatio(profile, opt.maxLag);
    [~, autocorrDecay, mean_autocorr, autocorrlm] = autocorrelationDecay(ST, lengthST, opt.T_ref);
end