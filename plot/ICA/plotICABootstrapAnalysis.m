function plotICABootstrapAnalysis(region, opt)
arguments
  region
  opt.significanceLevel = 0.05
  opt.dicoStep = 10
  opt.min_decade = 1.5
  opt.base = 1.2
  opt.threshold = 4
end
    [~, ~, ~, Ealpha, Ebeta, Egamma1, Egamma2, Egamma3] = ICABootstrapAnalysis(region, significanceLevel=opt.significanceLevel, dicoStep=opt.dicoStep, min_decade=opt.min_decade, base=opt.base, threshold=opt.threshold);
    % Combine all the variables into a matrix for easy access
    plotBoostrapAnalysis(Ealpha, Ebeta, Egamma1, Egamma2, Egamma3)
end