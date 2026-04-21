function p = kolmogorovSmirnovGoodnessFit(x,xmin,xmax,alpha)

    X = x((x>=xmin) & (x<=xmax));
    n = length(X);

    if n < 1000
        p = 0;

    else
        t = (xmin : xmax).';
        C = 1 / sum(t.^(-alpha));
        cdf_theoretical = cumsum(C*t.^(-alpha));
        [~,p] = kstest(X,'CDF',[t.',cdf_theoretical.']);

    end
end