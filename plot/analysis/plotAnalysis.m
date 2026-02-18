function plotAnalysis(T, alpha, Txmin, Txmax, Tp, S, beta, Sxmin, Sxmax, Sp, Alm, gam, x_shape, T_shape, shape, shape_gam, branching, crossCorr, branchingRatiolm, autocorrDecay, mean_autocorr, autocorrlm, T_ref)
    figure();
    subplot(2, 3, 1);
    plotDiscretePowerLawDensityFit(T, Txmin, Txmax, alpha);
    title(sprintf('lifetime PDF (p-value: {%.3f})', Tp));
    xlabel('lifetime');
    subplot(2, 3, 2);
    plotDiscretePowerLawDensityFit(S, Sxmin, Sxmax, beta, 'beta');
    title(sprintf('Size PDF (p-value: {%.3f})', Sp));
    xlabel('size');
    subplot(2, 3, 3);
    plotPowerLawFunction(Alm, expected_gamma=gam);
    title(sprintf('Area vs Lifetime fit (R² = {%.3f})', Alm.Rsquared.Ordinary));
    xlabel('lifetime');
    ylabel('Area');
    subplot(2, 3, 4);
    plotScalledShapeCollapseTransform(x_shape, T_shape, shape, shape_gam);
    subplot(2, 3, 5);
    plotBranchingRatioFit(branching, crossCorr, branchingRatiolm);
    subplot(2, 3, 6);
    plotAutocorrelationDecay(autocorrDecay, mean_autocorr, autocorrlm, T_ref);
end