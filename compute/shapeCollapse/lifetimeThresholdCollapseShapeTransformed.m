function [threshold_T, threshold_shape] = lifetimeThresholdCollapseShapeTransformed(T, shape, threshold)
    %keeps only avalanche with lifetime above the threshold
    indices = T >= threshold;
    threshold_T = T(indices);
    threshold_shape = shape(:, indices);
end

