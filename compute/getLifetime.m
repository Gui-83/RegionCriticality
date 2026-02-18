function T = getLifetime(aval_indeces)
    % compute the lifetime with aval_indeces (n, 2) which contains the start and end index of avalanches
    T = aval_indeces(:,2) - aval_indeces(:,1) + 1;
end

