function [S, T, ST, Ealpha, Ebeta, Egamma1, Egamma2, Egamma3] = rawRegionBootsrapAnalysis(regions, opt)
%wrapper for bootstrapAnalysis, it just extract data from region and feed
%it to the function. Return the data on which the boostrap was performed
%and the list of estimated alpha, gamma (for power law fit) and gammas (the
%3 different ways = scaling relation, shape collapse, Area vs Lifetime
%relation)
arguments
  regions
  opt.significanceLevel = 0.05
  opt.dicoStep = 10
  opt.min_decade = 1
  opt.base = 1.5
  opt.threshold = 4
  opt.samples = 0.01
  opt.bootstrap = 100
  opt.states = 'all'
  opt.regions = [-1]
end
    S = getRawRegionSize(regions, states=opt.states, regions=opt.regions);
    T = getRawRegionLifetime(regions, states=opt.states, regions=opt.regions);
    
    timeDependendentSize = getRawRegionSizeTimeDependent(regions, states=opt.states, regions=opt.regions);
    ST = separateAvalSizeTimeDependent(timeDependendentSize);

    [S, T, ST, Ealpha, Ebeta, Egamma1, Egamma2, Egamma3] = bootstrapAnalysis(S, T, ST, significanceLevel=opt.significanceLevel, dicoStep=opt.dicoStep, min_decade=opt.min_decade, base=opt.base, threshold=opt.threshold, samples=opt.samples, bootstrap=opt.bootstrap);
end