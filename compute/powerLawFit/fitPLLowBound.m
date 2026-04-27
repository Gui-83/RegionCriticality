function [alpha,xmin,xmax,p,h,D] = fitPLLowBound(x,n_search,significanceLevel,dicoStep)
% x: data
% significanceLevel: significance level for the Kolmogorov-Smirnov test
% dicoStep: number of step to perform the dichotomy

arguments
  x
  n_search (1,1) = 30
  significanceLevel (1,1) = 0.15
  dicoStep (1,1) = 100
end

x = nonUnique(x,3); % keep only data points which repeat at lest 3 times
[xminBound,xmaxBound] = bounds(x);

% build range for grid search
x_min = unique(ceil(logspace(log10(xminBound),log10(xmaxBound-1),n_search)));

% fit exponent, compute KS statistic
alpha = arrayfun(@(a) DiscreteBoundedPowerLawMLE(x,a,xmaxBound,dicoStep), x_min);
n_elem_min = 100;
D = arrayfun(@(a,b) ksStatistic(x,a,xmaxBound,b,n_elem_min), x_min, alpha);
[D,target] = min(D);

% reject range if not described by a power law
n_boot = 1000;
p = bootstrapKSp(x,x_min(target),xmaxBound,D,alpha(target),n_boot,dicoStep);
alpha = alpha(target);
xmin = x_min(target);
xmax = xmaxBound;
h = p > significanceLevel;

end


% --- helper functions ---


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


function [pmf,t] = powerLawPdf(x_min,x_max,alpha)
  t = (x_min : x_max).';
  W = exp(-alpha * log(t));
  pmf = W / sum(W);
end


function plotInDebug()
  % this function is not meant to be called, rather its code can be executed in debug mode to produce plots

  % plot empirical and fitted pdfs for a given range
  k = target;
  figure, plotDistr(x,'log',true)
  %XX = x(x >= x_min(k) & x <= xmaxBound); % plot pdf(x) only in range
  %figure, plotDistr(XX,'log',true)
  xline([x_min(k),xmaxBound],'k--')
  [pmf,t] = powerLawPdf(x_min(k),xmaxBound,alpha(k));
  pmf_data = histcounts(x,[t-0.5;t(end)+0.5],'Normalization','pdf').';
  pmf = pmf * trapz(t,pmf_data) / trapz(t,pmf); % rescale pdf to match whole data plot
  loglog(t,pmf)
  clearvars t pmf pmf_data

end