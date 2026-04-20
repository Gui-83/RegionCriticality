function y = preprocess(x)
    % remove outliers that appear only once
    uniqueVals = unique(x);
    [counts, ~] = histcounts(x, [uniqueVals; max(uniqueVals) + 1]);
    moreThanOnce = uniqueVals(counts > 1);
    y = x(ismember(x, moreThanOnce));
end