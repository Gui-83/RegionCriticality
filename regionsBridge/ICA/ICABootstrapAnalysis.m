function [S, T, ST, Ealpha, Ebeta, Egamma1, Egamma2, Egamma3] = ICABootstrapAnalysis(region, opt)
arguments
  region
  opt.significanceLevel = 0.05
  opt.dicoStep = 10
  opt.min_decade = 1
  opt.base = 1.5
  opt.threshold = 4
  opt.samples = 0.1
  opt.bootstrap = 1000
end
    S = catICAAvalSize(region);
    T = catICAAvalLifetime(region);
    
    timeDependendentSize = catICAAvalSizeTimeDependent(region);
    ST = separateAvalSizeTimeDependent(timeDependendentSize);

    [S, T, ST, Ealpha, Ebeta, Egamma1, Egamma2, Egamma3] = bootstrapAnalysis(S, T, ST, significanceLevel=opt.significanceLevel, dicoStep=opt.dicoStep, min_decade=opt.min_decade, base=opt.base, threshold=opt.threshold, samples=opt.samples, bootstrap=opt.bootstrap);
end