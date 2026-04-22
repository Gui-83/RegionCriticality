function lm = fitPowerFunction(y)

x = 1 : length(y);

non_zero = y > 0;
y = y(non_zero);
x = x(non_zero);

lm = fitlm(log(x),log(y));