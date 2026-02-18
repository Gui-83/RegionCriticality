function alpha = DiscreteBoundedPowerLawMLE(x, xmin, xmax, step)
    %perform MLE maximization to find the best power law parameter on the
    %raw data x and with the constrain xmin and xmax
    X = x((x>=xmin) & (x<=xmax));
    n = length(X);
    l1 = sum(log(X));
    t = xmin:xmax;
    log_t = log(t);

    upper = 6;
    lower = 1.1;
    
    %the dicotomie looks for the 0 of the log Likelihood derivative (see
    %docs/powerLawFit.md)
    for i = 1:step
        mid = (upper + lower)/2;
        D = DiscreteBoundedPowerLawLikelihoodDerivative(n, l1, t, log_t, mid);
        if D > 0
            lower = mid;
        else
            upper = mid;
        end
    end
    alpha = (upper + lower)/2;
end

