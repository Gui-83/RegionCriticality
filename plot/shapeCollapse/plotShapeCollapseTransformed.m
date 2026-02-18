function plotShapeCollapseTransformed(x, shape)
%mostly useless function, it plots the unscalled shapes of avalanches.
%Because there are not scalled it is bad
%x and shape are the output of transformCollapseShape (see README.md)
% x - (1, 100) list
% shape - (100, n) list where n is the number of unique avalanche Liftime
    hold on
    n = size(shape, 2);
    for i = 1:n
        plot(x, shape(:, i));
    end
end

