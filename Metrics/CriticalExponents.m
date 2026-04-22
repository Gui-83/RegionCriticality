function [alpha, beta, expected_gamma, gam_area, gam_shape, chi, diffgam] = CriticalExponents(session,region,states,events,window,opt)

arguments
  session (1,1) string
  region (1,1) string = "nr"
  states (:,1) string = "all"
  events (:,1) string = "all"
  window (1,1) = 0.05
  opt.labels (:,1) string = states+"_"+events
  opt.threshold (1,1) {mustBeNumeric,mustBeNonnegative} = 30
  opt.S (1,:) double = NaN
  opt.T (1,:) double = NaN
end

% remove "all" and events starting with "^"
reverse = cellfun(@(x) x(1)=='^',events);
events = string(cellfun(@(x) x([x(1)~='^',true(1,numel(x)-1)]),events,'UniformOutput',false));
unique_events = unique(events(events~="all"));
events = string(cellfun(@(y) y{end}, cellfun(@(x) strsplit(x,'/'),events,'UniformOutput',false),'UniformOutput',false));

R = regions(session,'states',unique(states),'events',unique_events);
R = R.computeAvalanches(window,1,opt.threshold);
window = R.avalWindow(region);

for s = 1 : length(states)

  statename = opt.labels(s);

  % 1. intervals, S, T

  if any(isnan(opt.S))
    S.(statename) = round(R.avalSizes(states(s),region) * window);
  else
    S.(statename) = opt.S;
  end

  if any(isnan(opt.T))
    if reverse(s)
      event_intervals = SubtractIntervals(R.eventIntervals,R.eventIntervals(events(s)));
    elseif events(s) == "all"
      event_intervals = [];
    else
      event_intervals = R.eventIntervals(events(s));
    end

    S.(statename) = round(R.avalSizes(states(s),region,'restriction',event_intervals) * window);
    intervals = R.avalIntervals(states(s),region,'restriction',event_intervals);
    T.(statename) = round((intervals(:,2) - intervals(:,1)) / window);
    
  else
    T.(statename) = opt.T;
  end

  % 2. power laws on S and T

  [beta.(statename),xmin_S.(statename),xmax_S.(statename),p_S.(statename),h_S.(statename)] = fitPLLowBound(S.(statename)); % swipeMLEDiscretePowerLawBounds(S.(statename));
  [alpha.(statename),xmin_T.(statename),xmax_T.(statename),p_T.(statename),h_T.(statename)] = fitPLLowBound(T.(statename)); % swipeMLEDiscretePowerLawBounds(T.(statename));
  expected_gamma.(statename) = (alpha.(statename) - 1) / (beta.(statename) - 1);

  % 3. Area

  A = getArea(T.(statename),S.(statename));
  lm = fitPowerFunction(A);
  gam_area.(statename) = lm.Coefficients.Estimate(2);

end

% TEMP
for s = 1 : length(states)

  
  

  %Shape Collapse
  aval_timeDependentSize = R.avalSizeOverTime(states(s),region);
  [ST, lengthST] = separateAvalSizeTimeDependent(aval_timeDependentSize);

  ST_AVG = AvalAverageSizeTimeDependent(ST, lengthST);
  [x, shape, T_cs] = transformCollapseShape(ST_AVG);

  [T_cs, shape] = lifetimeThresholdCollapseShapeTransformed(T_cs, shape, 4);
  gam_shape.(statename) = fitCollapseShape(T_cs, shape, 10);

  % Branching ratio

  %[profile, time] = R.avalProfiles(state,region);
  %[m, r, lm_branching] = branchingRatio(profile, 12);

  % Autocorrelation Decay

  [ST, lengthST] = separateAvalSizeTimeDependent(aval_timeDependentSize);
  [~, decay, mean_autocorr, lm_autocorr] = autocorrelationDecay(ST, lengthST, 1/3, min_uniqueST_lengths=5);
  chi.(statename) = -lm_autocorr.Coefficients.Estimate(2);

  diffgam.(statename) = gam_area.(statename)-expected_gamma.(statename);

end

end

% --- Extra code to plot examples in debug mode ---

function plotInDebug()
  % this function is not meant to be called, rather its code can be executed in debug mode to produce plots

  figure
  PlotIntervals(event_intervals,'alpha',.7)
  RasterPlot([intervals(:,1),.5*ones(size(intervals,1),1)])
  RasterPlot([intervals(:,2),.5*ones(size(intervals,1),1)])
  PlotIntervals(R.eventIntervals('slownr'),'color',[1,1,0])

  %%
  [~,axs] = makeFigure('',upper(region)+" "+opt.threshold+"% "+window+" s",[1,5],'size',[25,5.3]);
  fields = string(fieldnames(S));
  for i = 1 : length(fields)
    plotDiscretePowerLawDensityFit(S.(fields(i)),xmin_S.(fields(i)),xmax_S.(fields(i)),beta.(fields(i)),'\beta','ax',axs(i),'n_bins',200)
    legend(axs(i),'off')
    tit = fields(i) + " " + numel(S.(fields(i))) + " " + round(beta.(fields(i)),2);
    if h_S.(fields(i))
      tit = tit + " *";
    end
    title(axs(i),tit)
  end
  clearvars axs fields i

  %%
  [~,axs] = makeFigure('',upper(region)+" T "+opt.threshold+"% "+window+" s",[1,5],'size',[25,5.3]);
  fields = string(fieldnames(T));
  for i = 1 : length(fields)
    plotDiscretePowerLawDensityFit(T.(fields(i)),xmin_T.(fields(i)),xmax_T.(fields(i)),alpha.(fields(i)),'\alpha','ax',axs(i),'n_bins',200)
    legend(axs(i),'off')
    tit = fields(i) + " " + numel(T.(fields(i))) + " " + round(alpha.(fields(i)),2);
    if h_T.(fields(i))
      tit = tit + " *";
    end
    title(axs(i),tit)
  end
  clearvars axs fields i

end