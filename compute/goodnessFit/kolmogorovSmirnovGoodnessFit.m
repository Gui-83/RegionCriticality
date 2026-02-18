function p = kolmogorovSmirnovGoodnessFit (x, xmin, xmax, alpha)
    X = x((x>=xmin) & (x<=xmax));
    n = length(X);
    if n < 1000
        p = 0;
    else
        t = xmin:xmax;
        C = 1/sum(t.^(-alpha));
        cdf_theoretical = cumsum(C*t.^(-alpha));
        cdf_empirical = zeros(1, xmax-xmin+1);
        for i = xmin:xmax
            cdf_empirical(i-xmin+1) = sum((X <= i))/n;
        end
        d = max(abs(cdf_theoretical - cdf_empirical));
        p = 2*sum(arrayfun(@(k) (-1)^(k-1) * exp(-2 * k^2 * n * d^2), 1:10));
    end
end