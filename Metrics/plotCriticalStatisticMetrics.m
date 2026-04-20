function plotCriticalStatisticMetrics(region, state, S, T, xmin_S, xmin_T, xmax_S, xmax_T, alpha, beta, expected_gamma, lm, x, T_cs, shape, gam, m, r, lm_branching, decay, mean_autocorr, lm_autocorr)

%[fig,axs] = makeFigure('Critical Brain Statistics',sprintf('Critical Brain Statistics in the %s for the %s state', region, state),[2,4]);
[fig,axs] = makeFigure('Critical Brain Statistics',sprintf('Critical Brain Statistics in the %s for the %s state', region, state),[2,3],'format','poster');

plotDiscretePowerLawDensityFit(T, xmin_T, xmax_T, alpha, "alpha", ax=axs(1), title = "Power Law Distribution of T");
plotDiscretePowerLawDensityFit(S, xmin_S, xmax_S, beta, "beta", ax=axs(2), title = "Power Law Distribution of S");
plotPowerLawFunction(lm, ax=axs(3), expected_gamma=expected_gamma);
plotScalledShapeCollapseTransform(x, T_cs, shape, gam, ax=axs(4));
%plotBranchingRatioFit(m, r, lm_branching, ax=axs(5));
plotAutocorrelationDecay(decay, mean_autocorr, lm_autocorr, 1/3, ax = axs(5), loglog=false);
plotAutocorrelationDecay(decay, mean_autocorr, lm_autocorr, 1/3, ax = axs(6), loglog=true);

saveFig(gcf,"/mnt/hubel-data-103/Guillaume/NRAvalanche/Figures/svg_Files/Critical_Statistic_Metrics_ISR_(sans_branching_ratio)",'svg') % file_path: 'aaa/vvv/file_name'