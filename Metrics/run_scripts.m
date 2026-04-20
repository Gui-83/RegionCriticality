%session = "/mnt/hubel-data-131/perceval/Rat003_20231215/Rat003_20231215.xml";
session = "/mnt/hubel-data-131/perceval/Rat003_20231218/Rat003_20231218.xml";
R = regions(session,'states',["sws", "rem", "other"], 'events',"InfraSlowRhythm/slownr");
window = 0.05;

%R = R.computeAvalanches(window, 1, step=1);
R = R.computeAvalanches(window);

state = 'all';
region = "nr";

%plotAvalanchesCDF(R, region);

tmax=8000;
%nbz = 9;
%edges = linspace(0, tmax, nbz)';
%zones = [edges(1:end-1), edges(2:end)];
interval=[0, tmax];
st=80;
wd=2500;

testbin = 40;
slownr=0;
pt=0.05;

%plotCriticalityOverTime(R, region, "all", testbin, interval=interval, zones=zones);

[bin_centers, n_total, n_powerlaw,n_proportion, t_start, t_end, all_segments, seg_alpha, seg_beta, seg_gamma_exp, seg_gamma_area, seg_gamma_shape, seg_chi] = StateCriticalityOverTime_sliding_window(R, region, state, testbin, window, st=st, wd = wd, interval=interval, slownr=slownr, p_threshold=pt);
plotStateCriticalityOverTime_sliding_window(R, region, state, bin_centers, n_total, n_powerlaw,n_proportion, t_start, t_end, all_segments, seg_alpha, seg_beta, seg_gamma_exp, seg_gamma_area, seg_gamma_shape, seg_chi, wd, st, slownr=slownr);

%plotStateCriticalityOverTime_sliding_window(R, region, "all", testbin, interval=interval, st=st, wd=wd, slownr=0);

%plotStateCriticalityOverTime_sliding_window(R, region, "other", testbin, interval=interval, st=st, wd=wd, min_avals=50, p_threshold=0.0001);

%plotStateCriticalityOverTime_sliding_window(R, region, "sws", testbin, interval=interval, st=st, wd=wd, min_avals=10, p_threshold=0.005, slownr=1);

%plotStateCriticalityOverTime_sliding_window(R, region, "rem", testbin, interval=interval, st=st, wd=wd, min_avals=1000, p_threshold=0.05);

S = R.avalSizes(state,region) * window;

intervals = R.avalIntervals(state,region);

T = round((intervals(:,2) - intervals(:,1)) / window);

[beta, xmin_S, xmax_S, p_S] = swipeMLEDiscretePowerLawBounds(S, 0.05, 10, 0.1, 3, true);

%plotDiscretePowerLawDensityFit(S, xmin, xmax, beta, "beta");

[alpha, xmin_T, xmax_T, p_T] = swipeMLEDiscretePowerLawBounds(T, 0.05, 10, 1, 3, true);

%plotDiscretePowerLawDensityFit(T, xmin, xmax, alpha, "alpha");



%mask = (T >= 3) & (T <= 12);
%T = T(mask);
%S = S(mask);
A = getArea(T, S);

lm = fitPowerFunction(A);

%plotPowerLawFunction(lm);

aval_timeDependentSize = R.avalSizeOverTime(state,region);


[ST, lengthST] = separateAvalSizeTimeDependent(aval_timeDependentSize);

%mask = (lengthST >= 3) & (lengthST <= 12);
%ST = ST(mask);
%lengthST = lengthST(mask);

ST_AVG = AvalAverageSizeTimeDependent(ST, lengthST);
[x, shape, T_cs] = transformCollapseShape(ST_AVG);

[T_cs, shape] = lifetimeThresholdCollapseShapeTransformed(T_cs, shape, 4);
gam = fitCollapseShape(T_cs, shape, 10);

%plotScalledShapeCollapseTransform(x, T_cs, shape, gam);

[profile, time] = R.avalProfiles("sws",region);
[m, r, lm_branching] = branchingRatio(profile, 10);

%plotBranchingRatioFit(m, r, lm_branching);

[ST, lengthST] = separateAvalSizeTimeDependent(aval_timeDependentSize);

[~, decay, mean_autocorr, lm_autocorr] = autocorrelationDecay(ST, lengthST, 1/3, min_uniqueST_lengths=5);

%plotAutocorrelationDecay(decay, mean_autocorr, lm_autocorr, 1/3);

[fig,axs] = makeFigure('Critical Brain Statistics',sprintf('Critical Brain Statistics in the %s for the %s state', region, state),[2,3]);
plotDiscretePowerLawDensityFit(T, xmin_T, xmax_T, alpha, "alpha", ax=axs(1), title = "Power Law Distribution of T");
plotDiscretePowerLawDensityFit(S, xmin_S, xmax_S, beta, "beta", ax=axs(2), title = "Power Law Distribution of S");
plotPowerLawFunction(lm, ax=axs(3));
plotScalledShapeCollapseTransform(x, T_cs, shape, gam, ax=axs(4));
plotBranchingRatioFit(m, r, lm_branching, ax=axs(5));
plotAutocorrelationDecay(decay, mean_autocorr, lm_autocorr, 1/3, ax = axs(6));

