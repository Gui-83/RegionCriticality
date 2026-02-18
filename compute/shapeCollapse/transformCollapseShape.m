function [x, S, T] = transformCollapseShape(shape)
    %interpolate all the array in shape so that they all have the length m
    %x represent linear time between 0 and 1 for each shape
    n = numel(shape);
    m = 100;
    x = linspace(0, 1, m);
    S = zeros(m, n-1);
    T = zeros(1, n-1);
    for i = 1:n-1
        s = shape{i+1};
        t = length(s);
        z = linspace(0, t, t);
        s = interp1(z, s, x*t);
        S(:, i) = s;
        T(i) = t;
    end
end

