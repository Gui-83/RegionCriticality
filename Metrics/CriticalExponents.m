function [alpha, beta, expected_gamma, gam_area, gam_shape, chi, diffgam] = CriticalExponents(session, region, states, event, window, opt)

arguments
    session (1,1)
    region (1,1) string = "nr"
    states (1,:) string = ["all"]
    event (1,:) string = ["all"]
    %bin_size (1,1) {mustBeNumeric, mustBePositive} = 0.5
    window (1,1) = 0.05
    opt.S (1,:) double = [NaN, NaN]
    opt.T (1,:) double = [NaN, NaN]
end

% intervals, S, T

R = regions(session, 'states', ["sws", "rem", "other"],'events',"InfraSlowRhythm/slownr");
R = R.computeAvalanches(window);

for s = 1:length(states)

statename=states(s);

if any(isnan(opt.S))
    S = R.avalSizes(states(s),region) * window;
else
    S = opt.S;
end

if any(isnan(opt.T))
    intervals = R.avalIntervals(states(s),region);
    if event(s) == "slownr"
        intevent = R.eventIntervals(event(s));
        S = R.avalSizes(states(s),region, 'restriction',intevent) * window;
        intervals = R.avalIntervals(states(s),region, 'restriction',intevent);
        statename = "ISR";
    end
    if event(s) == "exceptslownr"
        intevent = SubtractIntervals(R.eventIntervals("all"), R.eventIntervals("slownr"));
        S = R.avalSizes(states(s),region, 'restriction',intevent) * window;
        intervals = R.avalIntervals(states(s),region, 'restriction',intevent);
        statename = "sws_non_ISR";
    end
    T = round((intervals(:,2) - intervals(:,1)) / window);
else
    T = opt.T;
end


%Power Laws on S and T

%display(statename);

[beta.(statename), xmin_S, xmax_S, p_S] = swipeMLEDiscretePowerLawBounds(S, 0.15, 10, 0.1, 3, true);

[alpha.(statename), xmin_T, xmax_T, p_T] = swipeMLEDiscretePowerLawBounds(T, 0.15, 10, 1, 3, true);

expected_gamma.(statename) = (alpha.(statename) - 1) / (beta.(statename) - 1);

%Area

A = getArea(T, S);
lm = fitPowerFunction(A);
gam_area.(statename) = lm.Coefficients.Estimate(2);
aval_timeDependentSize = R.avalSizeOverTime(states(s),region);

%Shape Collapse

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
