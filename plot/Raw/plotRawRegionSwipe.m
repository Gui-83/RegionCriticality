function plotRawRegionSwipe(regions, opt)
%wrapper for rawRegionSwipe it plots different criticality score for each
%state and region combination
arguments
  regions
  opt.significanceLevel = 0.05
  opt.dicoStep = 20
  opt.min_decade = 1
  opt.base = 2
  opt.threshold = 4
  opt.percentile = 60;
end
    [gammas, branchings, chis, decades] = rawRegionSwipe(regions, significanceLevel=opt.significanceLevel, dicoStep=opt.dicoStep, min_decade=opt.min_decade, base=opt.base, threshold=opt.threshold);
    v = var(gammas, 0, 3);
    precision = 1./v;
    h = figure();
    set(h, 'Name', 'Criticality score');
    subplot(2, 2, 1);
    plotCriticalityScore(regions.ids, regions.states, precision, logscale=true); %precision of the gammas estimation
    title('Gamma precision');
    subplot(2, 2, 2);
    plotCriticalityScore(regions.ids, regions.states, branchings); %branching ratio (if close to 1 then good)
    title('Branching ratio');
    subplot(2, 2, 3);
    decay = 1./chis;
    plotCriticalityScore(regions.ids, regions.states, decay); %slow decay of autocorelation (if big then good)
    title('Autocorrelation decay');
    subplot(2, 2, 4);
    plotCriticalityScore(regions.ids, regions.states, decades); %decade of which the power law were fited (if big then good)
    title('Size decade fit');
end