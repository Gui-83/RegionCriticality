function p = chiSquaredGodnessFit (x, xmin, xmax, alpha)
    X = x((x>=xmin) & (x<=xmax));
    n = length(X);
    if n < 1000
        p = 0;
    else
        t = xmin:xmax;
        C = 1/sum(t.^(-alpha));
        expected = C*t.^(-alpha);
        observed = zeros(1, xmax - xmin + 1);
        for i = xmin:xmax
            observed(i-xmin + 1) = sum(x == i)/n;
        end
        chi_square_stat = sum((observed - expected).^2 ./ expected);
        df = length(observed) - 1;
        p = 1 - chi2cdf(chi_square_stat, df);
    end
end