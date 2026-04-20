 function plotAutocorrelationDecay(autocorr, mean_autocorr, lm, T_ref, opt)
    % PLot autocorrelation decay (see Yang Tian Theoretical foundations of
    % studying)
    % autocorr : cell of autotoceralation function (the size of this cell
    % is equal to the number of unique Liftime). This variable represente
    % the left side of the equation 53 of the mentioned paper
    %criticality in the brain

    arguments
        autocorr
        mean_autocorr
        lm
        T_ref
        opt.ax = []
        opt.title = "Slow Decay of Auto-correlation"
        opt.loglog (1,1) logical = false
    end
    if isempty(opt.ax)
        opt.ax = gca;
    end


    n = length(autocorr);
    if opt.loglog
        set(opt.ax, 'XScale', 'log', 'YScale', 'log');
    end
    hold(opt.ax, 'on');

    %plots all the autocorr
    for i =1:n
        y = autocorr{i};
        y = y(y>0);
        if ~isempty(y)
            %x = linspace(0, 1, length(y));
            if opt.loglog
                x = linspace(1/length(y), 1, length(y));
                loglog(opt.ax, x, y, 'LineWidth', 1, 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
            else
                x = linspace(0, 1, length(y));
                semilogy(opt.ax, x, y, 'LineWidth', 1, 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
            end
        end
    end

    %plot the mean value of all autocorr
    %z = linspace(0, 1, length(mean_autocorr));
    if opt.loglog
        z = linspace(1/length(mean_autocorr), 1, length(mean_autocorr));
        mean_autocorr = mean_autocorr(mean_autocorr>0);
        loglog(opt.ax, z, mean_autocorr, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    else
        z = linspace(0, 1, length(mean_autocorr));
        semilogy(opt.ax, z, mean_autocorr, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end

    %plots the fit (equation 53 of the paper)
    chi = -lm.Coefficients.Estimate(2);
    r = lm.Coefficients.Estimate(1);

    fitindices = (z > T_ref);

    xfit = z(fitindices) - T_ref;
    yfit = xfit.^(-chi)*exp(r);

    if opt.loglog
        indices = z > T_ref;
        x_log = z(indices);
        y_log = mean_autocorr(indices);

        % filter non-positive values
        valid = y_log > 0;
        x_log = x_log(valid);
        y_log = y_log(valid);

        log_x = log(x_log);
        log_y = log(y_log);
        n_pts = length(log_x);

        best_r2 = -Inf;
        best_p  = [];
        best_range = [1, n_pts];

        % swipe over all [i,j] sub-intervals with minimum length
        min_pts = max(5, floor(n_pts * 0.2)); % au moins 20% des points
        for i = 1:n_pts
            for j = i+min_pts : n_pts
                hx = log_x(i:j);
                hy = log_y(i:j);
                p_fit = polyfit(hx, hy, 1);
                y_pred = polyval(p_fit, hx);
                SS_res = sum((hy - y_pred).^2);
                SS_tot = sum((hy - mean(hy)).^2);
                if SS_tot == 0, continue, end
                n = length(hx);
                r2_adj = 1 - (SS_res/SS_tot) * (n-1)/(n-2);
                if r2_adj > best_r2
                    best_r2   = r2_adj;
                    best_p    = p_fit;
                    best_range = [i, j];
                end
            end
        end

        % plot only the best sub-interval
        idx = best_range(1):best_range(2);
        chi_exp = -best_p(1);
        y_reg = exp(polyval(best_p, log_x(idx)));

        loglog(opt.ax, x_log, y_log, 'b.', 'HandleVisibility', 'off');
        loglog(opt.ax, x_log(idx), y_reg, 'r-', 'LineWidth', 2, ...
            'DisplayName', sprintf('\\chi_{exp} = %.3f (R²=%.2f)', chi_exp, best_r2));
        xlim(opt.ax, [T_ref, 1]);
    else
        semilogy(opt.ax, xfit + T_ref, yfit, 'LineWidth', 2, 'DisplayName', sprintf('\\chi = %.3f', chi));
    end

    if opt.loglog
        set(opt.ax, 'XScale', 'log', 'YScale', 'log');
    end
    
    ylo = min(mean_autocorr)/2;
    yhi = max(mean_autocorr)*2;
    if isfinite(ylo) && isfinite(yhi) && ylo < yhi
        ylim(opt.ax, [ylo, yhi]);
    end
    xline(opt.ax, T_ref, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('reference : %.2f', T_ref));

    xlabel(opt.ax, 't/T');
    ylabel(opt.ax, 'cross correlation')

    title(opt.ax, opt.title);
    legend(opt.ax, 'show');
    %plot(opt.ax,(1:100).^(-2));
end