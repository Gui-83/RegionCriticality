function scale_shape = scaleCollapseShape(shape, T, gam)
% scale the shape according to equation 49 of Yang Tian Theoretical foundations of
    % studying
    T = reshape(T, 1, []);
    scale_shape = shape.*T.^(1-gam);
end

