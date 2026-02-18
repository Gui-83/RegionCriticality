function plotRawRegionBoostrapAnalysis(regions,opt)
%wrapper for rawRegionBootsrapAnalysis
arguments
  regions
  opt.significanceLevel = 0.05
  opt.dicoStep = 10
  opt.min_decade = 1.5
  opt.base = 3
  opt.threshold = 4
  opt.states = 'all'
  opt.regions = [-1]
end
    [~, ~, ~, Ealpha, Ebeta, Egamma1, Egamma2, Egamma3] = rawRegionBootsrapAnalysis(regions, significanceLevel=opt.significanceLevel, dicoStep=opt.dicoStep, min_decade=opt.min_decade, base=opt.base, threshold=opt.threshold, states=opt.states, regions=opt.regions);
    % Combine all the variables into a matrix for easy access
    plotBoostrapAnalysis(Ealpha, Ebeta, Egamma1, Egamma2, Egamma3)
end

