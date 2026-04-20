function plotScalledShapeCollapseTransform(x, T, shape, gam, opt)
    %plot the shape collapse, if it works well all the shape in the shape
    %atgument should fall onto the same line (idealy a quadratic)
    % x - (1, 100) list
    % shape - (100, n) list where n is the number of unique avalanche Liftime
    % T - (1, n) list containing the corresponding lifetime of shape.

    arguments
        x
        T
        shape
        gam
        opt.ax = []
    end
    if isempty(opt.ax)
        opt.ax = gca;
    end
    
    % Calculate scaled and collapsed shape
    scale_shape = scaleCollapseShape(shape, T, gam);
    
    % Get number of columns (shapes) to plot
    n = size(scale_shape, 2);
    
    % Set up the plot
    hold(opt.ax, 'on');
    title(opt.ax, sprintf('Scaled Shape Collapse (we considered %d shapes (mean profiles))', n));
    xlabel(opt.ax, 't/T');
    ylabel(opt.ax, 'Scaled Shape');
    grid(opt.ax, 'on');
    
    % Plot each shape
    for i = 1:n
        plot(opt.ax, x, scale_shape(:,i), 'LineWidth', 1, 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off'); % Plot shapes in gray
    end
    
    % Compute the mean shape and plot it
    mean_shape = mean(scale_shape, 2);
    plot(opt.ax, x, mean_shape, 'LineWidth', 2, 'Color', 'b', 'DisplayName', sprintf('\\gamma = %.2f', gam)); % Plot mean shape in blue
    
    % Add a legend only for the mean shape
    legend(opt.ax, "show");
    
    % Ensure the plot is displayed correctly
    hold(opt.ax, 'off');
end