function [alpha, xmin, xmax, p] = swipeMLEDiscretePowerLawBounds(x, significanceLevel, dicoStep, min_decade, base, verbose)
    % x: data
    % significanceLevel: significance level for the Kolmogorov-Smirnov test
    % dicoStep: number of step to perform the dichotomy
    % min_decade: minimum number of decade between xmin and xmax
    % base: base of the grid search
    % verbose: display progress bar

    %see docs/powerLawFit.md for more information
    x = preprocess(x);
    xminBound = min(x);
    xmaxBound = max(x);

    if xmaxBound - xminBound - 10^min_decade < 0
        p = 0;
        alpha = NaN;
        xmin = 1;
        xmax = 1;
        return
    end
    
    current_p = 0;
    current_alpha = NaN;
    current_xmin = 1;
    current_xmax = 1;
    
    if verbose
        f = waitbar(0,'Please wait...');
    end

    total_iterations = get_total_iteration(xmaxBound, xminBound, min_decade, base);

    imax = upper_i(xmaxBound, xminBound, min_decade, base);
    iter = 0;
    for i = 0:imax
        xmin = getXmin(xminBound, base, i);
        [lowerj, upperj] = bound_j(xmaxBound, xmin, min_decade, base);
        for j = lowerj:upperj
            xmax = getXmax(xmin, base, j);
            iter = iter + 1;
            if verbose
                waitbar((iter-1)/total_iterations, f, ...
                    sprintf('Computing (%d/%d)', iter, total_iterations));
            end
            alpha = DiscreteBoundedPowerLawMLE(x, xmin, xmax, dicoStep); % perform MLE
            p = kolmogorovSmirnovGoodnessFit(x, xmin, xmax, alpha); % perform Kolmogorov-Smirnov test
            if (p > significanceLevel)
                if log10(xmax/xmin) > log10(current_xmax/current_xmin) % if the range is better that the previous fit, keeps the current fit
                    current_p = p;
                    current_alpha = alpha;
                    current_xmin = xmin;
                    current_xmax = xmax;
                end
            end
        end
    end
    alpha = current_alpha;
    xmin = current_xmin;
    xmax = current_xmax;
    p = current_p;
    if verbose
        close(f);
    end
end

function xmin = getXmin(xminBound, base, i)
    xmin = floor(xminBound + base^i);
end

function xmax = getXmax(xmin, base, j)
    xmax = floor(xmin + base^j);
end

function i = upper_i(xmaxBound, xminBound, min_decade, base)
    %see docs/powerLawFit.md for more information
    i = floor(log(xmaxBound - xminBound - 10^min_decade)/log(base));
end

function [lower, upper] = bound_j(xmaxBound, xmin, min_decade, base)
    %see docs/powerLawFit.md for more information
    lower = floor(log(10^min_decade)/log(base));
    upper = floor(log(xmaxBound - xmin)/log(base));
end

function total = get_total_iteration(xmaxBound, xminBound, min_decade, base)
    imax = upper_i(xmaxBound, xminBound, min_decade, base);
    total = 0;
    for i = 0:imax
        xmin = getXmin(xminBound, base, i);
        [lowerj, upperj] = bound_j(xmaxBound, xmin, min_decade, base);
        total = total + upperj - lowerj +1;
    end
end
    
