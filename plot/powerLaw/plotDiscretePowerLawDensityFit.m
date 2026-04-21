function plotDiscretePowerLawDensityFit(x,xmin,xmax,alpha,fitVariable,opt)
%plot empirical powerlaw density of x and the fited bounded powerlaw
% x - 1-D list raw data
% xmin and xmax : bound of the fit powerlaw
% alpha - parameter of the power law
%fitVariable - either alpha or beta (name of the variable displayed on
%the plot
arguments
  x
  xmin
  xmax
  alpha
  fitVariable (1,1) string = "\alpha"
  opt.ax = gca
  opt.title = ""
end

x = nonUnique(x);
plotDistr(x,'log',true,'label','empirical pdf','ax',opt.ax)

t = (xmin : xmax).';
w = exp(-alpha * log(t));
C = 1 / sum(w);
pmf = C * w;
loglog(opt.ax,t,pmf,'--','DisplayName',"fit: "+round(alpha,2)+", range: "+round(log10(xmax/xmin),2))

xline(opt.ax,[xmin,xmax],'k--','HandleVisibility', 'off')

% count avalanches between xmin and xmax
%n_in_range = sum(x >= xmin & x <= xmax);
%n_total_avals = length(x);

title(opt.ax, opt.title);
end