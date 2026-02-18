function plotICARegionSwipe(regions, opt)
arguments
  regions
  opt.significanceLevel = 0.2
  opt.dicoStep = 10
  opt.min_decade = 1
  opt.base = 1.5
  opt.threshold = 4
end
    gammas = ICARegionSwipe(regions, significanceLevel=opt.significanceLevel, dicoStep=opt.dicoStep, min_decade=opt.min_decade, base=opt.base, threshold=opt.threshold);
    v = var(gammas, 0, 3);
    figure;
    plotCriticalityScore(regions.ids(regions.ids > 0), regions.states, -log(v));
end