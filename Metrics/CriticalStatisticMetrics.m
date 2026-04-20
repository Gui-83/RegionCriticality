function [S, T, xmin_S, xmin_T, xmax_S, xmax_T, alpha, beta, expected_gamma, lm, x, T_cs, shape, gam, m, r, lm_branching, decay, mean_autocorr, lm_autocorr] = CriticalStatisticMetrics(R, region, state, event, window)

arguments
    R (1,1) regions
    region (1,1) string = "nr"
    state (1,1) string = "all"
    event (1,1) string = "all"
    %bin_size (1,1) {mustBeNumeric, mustBePositive} = 0.5
    window (1,1) = 0.05
end

% intervals, S, T

S = R.avalSizes(state,region) * window;

intervals = R.avalIntervals(state,region);
if event ~= "all"
        intevent = R.eventIntervals(event);
        S = R.avalSizes(state,region, 'restriction',intevent) * window;
        intervals = R.avalIntervals(state,region, 'restriction',intevent);
end
T = round((intervals(:,2) - intervals(:,1)) / window);


%Power Laws on S and T

[beta, xmin_S, xmax_S, p_S] = swipeMLEDiscretePowerLawBounds(S, 0.15, 10, 0.1, 3, true);

[alpha, xmin_T, xmax_T, p_T] = swipeMLEDiscretePowerLawBounds(T, 0.15, 10, 1, 3, true);

expected_gamma = (alpha-1)/(beta-1);

%Area

A = getArea(T, S);
lm = fitPowerFunction(A);
aval_timeDependentSize = R.avalSizeOverTime(state,region);

%Shape Collapse

[ST, lengthST] = separateAvalSizeTimeDependent(aval_timeDependentSize);

ST_AVG = AvalAverageSizeTimeDependent(ST, lengthST);
[x, shape, T_cs] = transformCollapseShape(ST_AVG);

[T_cs, shape] = lifetimeThresholdCollapseShapeTransformed(T_cs, shape, 4);
gam = fitCollapseShape(T_cs, shape, 10);

% Branching ratio

[profile, time] = R.avalProfiles(state,region);
[m, r, lm_branching] = branchingRatio(profile, 12);

% Autocorrelation Decay

[ST, lengthST] = separateAvalSizeTimeDependent(aval_timeDependentSize);
[~, decay, mean_autocorr, lm_autocorr] = autocorrelationDecay(ST, lengthST, 1/3, min_uniqueST_lengths=5);
