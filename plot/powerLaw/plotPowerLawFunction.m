function plotPowerLawFunction(lm, opt)
%This is used to plot the Area vs Lifetime relation. It takes a fitlm
%object as argument and pplot the regression on loglog
%if expected_gamma is specified it will be displayed on the plots. (the
%expected gamma iss given by the scaling relation)
arguments
    lm
    opt.expected_gamma = -1
end
    logx = lm.Variables{:,1};
    logy = lm.Variables{:,2};

    x = exp(logx);
    y = exp(logy);
    
    gam = lm.Coefficients.Estimate(2);
    intercept = lm.Coefficients.Estimate(1);
    a = exp(intercept);
    y_fit = a * x.^gam;
    
    
    scatter(x, y, 'filled', 'DisplayName', 'Original Data', 'MarkerFaceColor', 'blue');
    hold on;
    if opt.expected_gamma > 0
        leg = sprintf('\\gamma = {%.2f} (expected : \\gamma = {%.2f})', gam, opt.expected_gamma);
    else
        leg = sprintf('\\gamma = {%.2f}', gam);
    end
    plot(x, y_fit, '-r', 'LineWidth', 2, 'DisplayName', leg);
    
    set(gca, 'XScale', 'log', 'YScale', 'log');
    legend('show');
end
