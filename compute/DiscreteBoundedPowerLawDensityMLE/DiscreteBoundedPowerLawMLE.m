function alpha = DiscreteBoundedPowerLawMLE(x,xmin,xmax,step)
% perform MLE maximization to find the best power law parameter on data 'x', with the constraint xmin and xmax

X = x(x >= xmin & x <= xmax);
n = length(X);
l1 = sum(log(X));

t = (xmin : xmax).';
log_t = log(t);

lower = 1.1;
upper = 6;

% dichotomy to find 0 of the log-Likelihood derivative (see docs/powerLawFit.md)
for i = 1 : step
  mid = (upper + lower) / 2;

  % derivative
  log_w = -mid * log_t;
  w = exp(log_w - max(log_w)); % avoid overflow in exp
  l2 = sum(log_t .* w);
  l3 = sum(w);
  D = -l1 + n * l2 / l3;

  if D > 0
    lower = mid;
  else
    upper = mid;
  end
end

alpha = (upper + lower) / 2;