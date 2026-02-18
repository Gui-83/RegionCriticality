function lm = fitPowerFunction(y)
    n = length(y);
    x = 1:n;

    non_zero = y > 0;
    y = y(non_zero);
    x = x(non_zero);

    lm = fitlm(log(x), log(y));
end

