function [alpha,xmin,xmax,p] = swipeMLEDiscretePowerLawBounds(x,significanceLevel,dicoStep,min_decade,base)
% x: data
% significanceLevel: significance level for the Kolmogorov-Smirnov test
% dicoStep: number of step to perform the dichotomy
% min_decade: minimum number of decade between xmin and xmax
% base: base of the logarithmic grid search

arguments
  x
  significanceLevel (1,1) = 0.15
  dicoStep (1,1) = 100
  min_decade (1,1) = 0.1
  base (1,1) = 1.5
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

% build ranges for grid search, i: x min, j: x max, see docs/powerLawFit.md for info
i_min = double(base < 2);
i_max = floor(log(xmaxBound - xminBound - 10^min_decade) / log(base));
x_min = logCeil(xminBound,base,(i_min:i_max).',Inf);
j_min = floor(log(10^min_decade) / log(base));
j_min = j_min + double(j_min < 1 && base < 2);
j_max = ceil(log(xmaxBound - x_min) / log(base));
x_max = arrayfun(@(x,y) logCeil(x,base,j_min:y,xmaxBound), x_min, j_max, 'UniformOutput', false); % xmaxBound is always included in ranges
ranges = [repelem(x_min,cellfun(@numel,x_max)),[x_max{:}].'];

% fit exponent, copmute KS statistic
n_elem_min = 100;
alpha = arrayfun(@(a,b) DiscreteBoundedPowerLawMLE(x,a,b,dicoStep), ranges(:,1), ranges(:,2));
D = arrayfun(@(a,b,c) ksStatistic(x,a,b,c,n_elem_min), ranges(:,1), ranges(:,2), alpha);
valid = ~isnan(D);
ranges = ranges(valid,:);
D = D(valid);
alpha = alpha(valid);

% reject ranges not described by a power law
n_boot = 500;
p = arrayfun(@(a,b,c,d) bootstrapKSp(x,a,b,c,d,n_boot,dicoStep), ranges(:,1), ranges(:,2), D, alpha);
ok = find(p > significanceLevel);
[~,final] = max(diff(ranges(ok,:),1,2));

alpha = alpha(ok(final));
xmin = ranges(ok(final),1);
xmax = ranges(ok(final),2);
p = p(ok(final));

end


% --- helper functions ---


function y = logCeil(x,base,i,cap)
  y = floor(x + base.^i);
  y(y>cap) = cap;
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


function p = bootstrapKSp(x,xmin,xmax,D,alpha,n_boot,dicoStep)
  X = x(x >= xmin & x <= xmax);
  n = numel(X);

  t = (xmin : xmax).';
  w = exp(-alpha * log(t));
  C = 1 / sum(w);
  pmf = C * w;
  %cdf_theoretical = cumsum(pmf);

  D_boot = zeros(n_boot,1);
  for b = 1 : n_boot
    % sample from discrete power law via inverse CDF
    X_boot = randsample(t,n,true,pmf);
    % refit alpha
    alpha_boot = DiscreteBoundedPowerLawMLE(X_boot,xmin,xmax,dicoStep);
    % theoretical CDF
    C = 1 / sum(t.^(-alpha_boot));
    cdf_th = cumsum(C * t.^(-alpha_boot));
    % empirical CDF via counts
    cdf_emp = histcounts(X_boot,[t;xmax+1],'Normalization','cdf').';
    % bootstrapped statistic
    D_boot(b) = max(abs(cdf_emp - cdf_th));
  end

  p = mean(D_boot >= D);
end


% --- Extra code to plot examples in debug mode ---


function powerLawPdf()

end


function plotInDebug()
  % this function is not meant to be called, rather its code can be executed in debug mode to produce plots

  % plot empirical and fitted pdfs for a given range
  k = 10;
  XX = x(x >= ranges(k,1) & x <= ranges(k,2));
  figure, plotDistr(XX,'log',true)
  T = (ranges(k,1) : ranges(k,2)).';
  W = exp(-alpha(k) * log(T));
  C = 1 / sum(W);
  pmf = C * W;
  loglog(T,pmf)
  clearvars XX T W C pmf
  % I CAN FIND A COEFF to do pmf/coeff AND GET NICE PLOT

end