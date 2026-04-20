function alpha = DiscreteBoundedPowerLawMLE(x, xmin, xmax, step)
    % perform MLE maximization to find the best power law parameter on raw data x, with the constrain xmin and xmax
    X = x((x>=xmin) & (x<=xmax));
    n = length(X);
    l1 = sum(log(X));
    t = xmin:xmax;
    log_t = log(t);

    upper = 6;
    lower = 1.1;
    
    % the dichotomy looks for the 0 of the log Likelihood derivative (see docs/powerLawFit.md)
    for i = 1:step
        mid = (upper + lower)/2;
        D = derivative(n, l1, t, log_t, mid);
        if D > 0
            lower = mid;
        else
            upper = mid;
        end
    end
    alpha = (upper + lower)/2;
end

% --- helper functions ---

function D = derivative(n, l1, t, log_t, alpha)
    t_alpha = t.^(-alpha);
    l2 = sum(log_t.*t_alpha);
    l3 = sum(t_alpha);

    D = -l1 + n * l2 / l3;
end