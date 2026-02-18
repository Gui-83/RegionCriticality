function area = getArea(T, S)
    % compute the Area given the lifetime and size
    if isempty(T)
        area = [];
        return
    end
    [uniqueLifetime, ~, idx] = unique(T);
    meanSizes = accumarray(idx, S, [], @mean);
    area = zeros(max(T), 1);
    area(uniqueLifetime) = meanSizes;
end

