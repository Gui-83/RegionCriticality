function [alpha,xmin,xmax,p] = swipeMLEDiscretePowerLawBounds(x,significanceLevel,dicoStep,min_decade,base)
% x: data
% significanceLevel: significance level for the Kolmogorov-Smirnov test, e.g., 0.1
% dicoStep: number of step to perform the dichotomy
% min_decade: minimum number of decade between xmin and xmax
% base: base of the logarithmic grid search

arguments
  x
  significanceLevel (1,1) = 0.15
  dicoStep (1,1) = 10
  min_decade (1,1) = 0.1
  base (1,1) = 2
end

% see docs/powerLawFit.md for more information
x = nonUnique(x);
xminBound = min(x);
xmaxBound = max(x);

if xmaxBound - xminBound - 10^min_decade < 0
  p = 0;
  alpha = NaN;
  xmin = 1;
  xmax = 1;
  return
end

% build ranges for grid search, i: x min, j: x max
i_max = floor(log(xmaxBound - xminBound - 10^min_decade) / log(base)); % see docs/powerLawFit.md for info
x_min = logCeil(xminBound,base,(0:i_max).');
[lower_j,upper_j] = bound_j(xmaxBound,x_min,min_decade,base);
x_max = arrayfun(@(x,y,z) logCeil(x,base,y:z),x_min,lower_j,upper_j,'UniformOutput',false);
ranges = [repelem(x_min,cellfun(@numel,x_max)),[x_max{:}].'];

% fit exponent, copmute KS statistic
alpha = arrayfun(@(a,b) DiscreteBoundedPowerLawMLE(x,a,b,dicoStep), ranges(:,1), ranges(:,2));
D = arrayfun(@(a,b,c) ksStatistic(x,a,b,c,100), ranges(:,1), ranges(:,2), alpha);
valid = ~isnan(D);
ranges = ranges(valid,:);
D = D(valid);
alpha = alpha(valid);

% reject ranges not described by a power law
p = arrayfun(@(a,b,c,d) bootstrapKSp(x,a,b,c,d,1000), ranges(:,1), ranges(:,2), D, alpha);
ok = find(p > significanceLevel);
[~,final] = max(diff(ranges(ok,:),1,2));

alpha = alpha(ok(final));
xmin = ranges(ok(final),1);
xmax = ranges(ok(final),2);
p = p(ok(final));

end


% --- helper functions ---


function y = logCeil(x,base,i)
  y = floor(x + base.^i);
end


function [lower,upper] = bound_j(xmaxBound,xmin,min_decade,base)
  % see docs/powerLawFit.md for more information
  lower = repmat(floor(log(10^min_decade) / log(base)), size(xmin));
  upper = floor(log(xmaxBound - xmin) / log(base));
end


function D = ksStatistic(x,xmin,xmax,alpha,n_min)
  X = x(x >= xmin & x <= xmax);
  if length(X) < n_min
    D = NaN;
  else
    t = (xmin : xmax).';
    w = exp(-alpha * log(t));
    C = 1 / sum(w);
    pmf = C * w;
    cdf_theoretical = cumsum(pmf);
    cdf_emp = histcounts(X,[t;xmax+1],'Normalization','cdf').';
    D = max(abs(cdf_emp - cdf_theoretical));
  end
end


function p = bootstrapKSp(x,xmin,xmax,D,alpha,n_boot)
  X = x(x >= xmin & x <= xmax);
  n = numel(X);

  t = (xmin : xmax).';
  w = exp(-alpha * log(t));
  C = 1 / sum(w);
  pmf = C * w;
  cdf_theoretical = cumsum(pmf);

  D_boot = zeros(n_boot,1);
  for b = 1 : n_boot
    % sample from discrete power law via inverse CDF
    X_boot = randsample(t,n,true,pmf);
    % empirical CDF via counts
    cdf_emp = histcounts(X_boot,[t;xmax+1],'Normalization','cdf').';
    % bootstrapped statistic
    D_boot(b) = max(abs(cdf_emp - cdf_theoretical));
  end

  p = mean(D_boot >= D);
end