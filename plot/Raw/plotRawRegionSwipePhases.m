function plotRawRegionSwipePhases(session, opt)
%wrapper for rawRegionSwipePhases and plotCriticalityScorePhases. It plots
%the criticality score for all state, regions and phasses combination
arguments
  session
  opt.significanceLevel = 0.05
  opt.dicoStep = 10
  opt.min_decade = 1
  opt.base = 3
  opt.threshold = 4
end
    [Ids, States, Phases, G] = rawRegionSwipePhases(session, significanceLevel=opt.significanceLevel, dicoStep=opt.dicoStep, min_decade=opt.min_decade, base=opt.base, threshold=opt.threshold);
    plotCriticalityScorePhases(Ids, States, Phases, G);
end