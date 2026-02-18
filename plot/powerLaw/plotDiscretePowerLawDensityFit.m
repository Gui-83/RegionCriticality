function plotDiscretePowerLawDensityFit(x, xmin, xmax, alpha, fitVariable)
    %plot empirical powerlaw density of x and the fited bounded powerlaw
    % x - 1-D list raw data
    % xmin and xmax : bound of the fit powerlaw
    % alpha - parameter of the power law
    %fitVariable - either alpha or beta (name of the variable displayed on
    %the plot
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

    loglog(bin_centers, counts, 'bo', 'DisplayName', 'Empirical PDF');
    hold on;
    
    x_vals = xmin:xmax;
    C2 = 1/sum((xmin:xmax).^(-alpha));
    power_law_pdf = C2*x_vals.^(-alpha);
    loglog(x_vals, power_law_pdf, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Estimated Power Law PDF (%s = %.2f)', variableName, alpha));
    
    ymin = min(power_law_pdf)/10;
    ymax = max(power_law_pdf)*10;
    loglog([xmin, xmin], [ymin, ymax], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');  % Vertical line at xmin
    loglog([xmax, xmax], [ymin, ymax], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');  % Vertical line at xmax

    legend show;
end
