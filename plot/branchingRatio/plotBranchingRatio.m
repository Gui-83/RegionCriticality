function plotBranchingRatio(profile, kmax)
%wrapper for plotBranchingRatioFit. it first compute the branching from raw
%data then plot it
    [m, r, lm] = branchingRatio(profile, kmax);
    plotBranchingRatioFit(m, r, lm);
end

