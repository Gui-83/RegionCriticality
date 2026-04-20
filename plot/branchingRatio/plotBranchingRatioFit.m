function plotBranchingRatioFit(m, r, lm, opt)
    %plot the linear regression used to compute the branching ratio : Inferring collective
    %dynamical states from widely, Jens Wilting1 & Viola Priesemann
    %unobserved systems
    % m : is the estimated branching ratio
    % r : list of laged branching ratio
    % lm : fit of r against a power function
    arguments
        m
        r
        lm
        opt.ax = []
        opt.title = "Branching ratio"
    end
    if isempty(opt.ax)
        opt.ax = gca;
    end

    b = exp(lm.Coefficients.Estimate(1));
    kmax = length(r);
    x = 1:kmax;
    y = b*m.^x;
    semilogy(opt.ax, x, r, 'bo', 'HandleVisibility', 'off');
    hold(opt.ax, 'on');
    semilogy(opt.ax, x, y, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Estimated branching ratio : r = %.5f', m));
    xlabel(opt.ax, 'lags')
    ylabel(opt.ax, 'branching ratio')

    title(opt.ax, opt.title);
    set(opt.ax, 'YScale', 'log');
    legend(opt.ax, 'show');
end

