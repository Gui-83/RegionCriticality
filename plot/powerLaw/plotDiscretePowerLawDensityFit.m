function plotDiscretePowerLawDensityFit(x, xmin, xmax, alpha, fitVariable, opt)
    %plot empirical powerlaw density of x and the fited bounded powerlaw
    % x - 1-D list raw data
    % xmin and xmax : bound of the fit powerlaw
    % alpha - parameter of the power law
    %fitVariable - either alpha or beta (name of the variable displayed on
    %the plot
    arguments
        x
        xmin
        xmax
        alpha
        fitVariable (1,1) string = "alpha"
        opt.ax = []
        opt.title = "Power Law Distribution"
    end
    if isempty(opt.ax)
        opt.ax = gca;
    end
    x = preprocess(x);
    if nargin < 5
        fitVariable = 'alpha';
    end

    if strcmp(fitVariable, 'alpha')
        variableName = '\alpha';
    elseif strcmp(fitVariable, 'beta')
        variableName = '\beta';
    else
        error('fitVariable must be either "alpha" or "beta"');
    end
    
    n = max(x);
    counts = zeros(1, n);
    bin_centers = 1:n;
    for i = 1:n
        counts(i) = sum(x == i);
    end
    C1 = 1/sum(counts(1,xmin:xmax)); %normalisation constant
    indices = counts > 0;
    counts = counts(indices);
    bin_centers = bin_centers(indices);
    if (xmax - xmin > 0)
        counts = C1*counts;
    end

    loglog(opt.ax, bin_centers, counts, 'bo', 'DisplayName', 'Empirical PDF');
    hold(opt.ax, 'on');
    set(opt.ax, 'XScale', 'log', 'YScale', 'log');
    
    x_vals = xmin:xmax;
    C2 = 1/sum((xmin:xmax).^(-alpha));
    power_law_pdf = C2*x_vals.^(-alpha);
    loglog(opt.ax, x_vals, power_law_pdf, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Estimated Power Law PDF (%s = %.2f)', variableName, alpha));
    
    ymin = min(power_law_pdf)/10;
    ymax = max(power_law_pdf)*10;
    loglog(opt.ax, [xmin, xmin], [ymin, ymax], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');  % Vertical line at xmin
    loglog(opt.ax, [xmax, xmax], [ymin, ymax], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');  % Vertical line at xmax

% labels xmin / xmax sous les lignes
    text(opt.ax, xmin, ymin, sprintf('x_{min}=%g', xmin), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12);
    text(opt.ax, xmax, ymin, sprintf('x_{max}=%g', xmax), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 12);

    annotation_y = ymin * 1.2;
    plot(opt.ax, [xmin, xmax], [annotation_y, annotation_y], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(opt.ax, xmin*1.05, annotation_y, '<k', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(opt.ax, xmax*0.95, annotation_y, '>k', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'HandleVisibility', 'off');

% valeur log(xmax/xmin) au centre
    text(opt.ax, sqrt(xmin*xmax), ymin*1.2, sprintf('log range : %.2f', log10(xmax/xmin)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);

    % count avalanches between xmin and xmax
    n_in_range = sum(x >= xmin & x <= xmax);
    n_total_avals = length(x);

    plot(opt.ax, NaN, NaN, 'HandleVisibility', 'on', ...
    'DisplayName', sprintf("Proportion d'avalanches = %d/%d", n_in_range, n_total_avals), ...
    'LineStyle', 'none', 'Marker', 'none');
    
    title(opt.ax, opt.title);
    legend(opt.ax,'show');
end
