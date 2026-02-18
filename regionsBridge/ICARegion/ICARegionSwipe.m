function g = ICARegionSwipe(regions, opt)
arguments
  regions
  opt.significanceLevel = 0.2
  opt.dicoStep = 10
  opt.min_decade = 1
  opt.base = 1.5
  opt.threshold = 4
end
    g = regionsSwipeAnalysis(regions, @getICARegionSize, @getICARegionLifetime, @getICARegionArea, @getICARegionSizeTimeDependent, significanceLevel=opt.significanceLevel, dicoStep=opt.dicoStep, min_decade=opt.min_decade, base=opt.base, threshold=opt.threshold);
end