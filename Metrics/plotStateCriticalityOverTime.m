
function plotStateCriticalityOverTime(R, region, state, bin_size, window, opt)
arguments
    R (1,1) regions
    region (1,1) string
    state (1,1) string = "all"
    bin_size (1,1) {mustBeNumeric, mustBePositive} = 0.5
    window (1,1) = 0.05
    opt.interval (1,2) double = [NaN, NaN]
    opt.states (1,:) string = ["sws", "rem", "other"]
    opt.colors (:,3) double = [0.2 0.4 0.8; 0.8 0.2 0.2; 0.4 0.8 0.2]
    opt.min_avals (1,1) {mustBeNumeric} = 50
    opt.p_threshold (1,1) double = 0.05
    opt.n_swipe (1,1) double = 10
    opt.zones (:,2) double = []
end

% --- 1. Get all avalanches ---
intervals_all = R.avalIntervals('all', region);
S_all = R.avalSizes('all', region) * window;
T_all = round((intervals_all(:,2) - intervals_all(:,1)) / window);

% filter by state if requested
if state ~= "all"
    [~,~,s_index,~] = R.arrayInd(state, region);
    state_times = R.state.times{s_index};
    
    % filtrer
    [~, ind1] = Restrict(intervals_all(:,1), state_times);
    [~, ind2] = Restrict(intervals_all(:,2), state_times);
    ind = intersect(ind1, ind2);
    intervals_all = intervals_all(ind, :);
    S_all = S_all(ind);
    T_all = T_all(ind);
    
    % conserver les durées
    durations = intervals_all(:,2) - intervals_all(:,1);
    
    % shifter les starts
    aval_starts = Restrict(intervals_all(:,1), state_times, shift=true);
    
    % reconstruire les stops à partir des starts shiftés + durées
    aval_stops = aval_starts + durations;
else
    aval_starts = intervals_all(:,1);
    aval_stops  = intervals_all(:,2);
end

% --- 2. Build time bins ---
if any(isnan(opt.interval))
    t_start = min(aval_starts);
    t_end   = max(aval_stops);
else
    t_start = opt.interval(1);
    t_end   = opt.interval(2);
end

bin_edges = t_start : bin_size : t_end;

if bin_edges(end) < t_end
    bin_edges(end+1) = t_end;
end
n_bins = length(bin_edges) - 1;
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

% --- 3. Build segments for power law fitting ---
if ~isempty(opt.zones)
    all_segments = [opt.zones, zeros(size(opt.zones,1),1)];
else
    all_segments = [t_start, t_end, 0];
end
all_segments = sortrows(all_segments, 1);

% --- 4. Fit power law per segment ---
n_seg = size(all_segments, 1);

seg_alpha = NaN(n_seg, 1);
seg_beta = NaN(n_seg, 1);
seg_gamma_exp = NaN(n_seg, 1);
seg_gamma_area = NaN(n_seg, 1);
seg_gamma_shape = NaN(n_seg, 1);
seg_chi = NaN(n_seg, 1);

seg_xminS = NaN(n_seg, 1);
seg_xmaxS = NaN(n_seg, 1);
seg_xminT = NaN(n_seg, 1);
seg_xmaxT = NaN(n_seg, 1);

% Pre-compute outside the loop
aval_tDS_raw = R.avalSizeOverTime('all', region);
% aval_tDS_raw est un vecteur concaténé avec des 0 séparateurs
% il faut extraire les profils de toutes les avalanches d'abord
[ST_all, lenST_all] = separateAvalSizeTimeDependent(aval_tDS_raw);
% maintenant ST_all{i} correspond à l'avalanche i dans intervals_all original

% filter by state
if state ~= "all"
    ST_all    = ST_all(ind);    % ind = intersect(ind1,ind2) calculé plus tôt
    lenST_all = lenST_all(ind);
end
% maintenant length(ST_all) == length(aval_starts) par construction


for seg = 1:n_seg
    t1 = all_segments(seg, 1);
    t2 = all_segments(seg, 2);
    in_seg = aval_starts >= t1 & aval_stops <= t2;
    S_seg = S_all(in_seg);
    T_seg = T_all(in_seg);
    if sum(in_seg) < opt.min_avals
        %display trop peu;
        continue
    end

    % filter ST by segment
    ST_seg    = ST_all(in_seg);
    lenST_seg = lenST_all(in_seg);

    try
        [beta_seg, xmin_S, xmax_S] = swipeMLEDiscretePowerLawBounds(S_seg, opt.p_threshold, opt.n_swipe, 0.1, 3, false);
        [alpha_seg, xmin_T, xmax_T] = swipeMLEDiscretePowerLawBounds(T_seg, opt.p_threshold, opt.n_swipe, 1, 3, false);
        seg_xminS(seg) = xmin_S;
        seg_xmaxS(seg) = xmax_S;
        seg_xminT(seg) = xmin_T;
        seg_xmaxT(seg) = xmax_T;

        seg_alpha(seg) = alpha_seg;
        seg_beta(seg)  = beta_seg;
        seg_gamma_exp(seg) = (alpha_seg - 1) / (beta_seg - 1);

        % gamma_area
        try
            A_seg = getArea(T_seg, S_seg);
            lm_area = fitPowerFunction(A_seg);
            seg_gamma_area(seg) = lm_area.Coefficients.Estimate(2);
        catch, end

         % gamma_shape
        try
            ST_AVG_seg = AvalAverageSizeTimeDependent(ST_seg, lenST_seg);
            [x_seg, shape_seg, T_cs_seg] = transformCollapseShape(ST_AVG_seg);
            [T_cs_seg, shape_seg] = lifetimeThresholdCollapseShapeTransformed(T_cs_seg, shape_seg, 4);
            seg_gamma_shape(seg) = fitCollapseShape(T_cs_seg, shape_seg, 10);
        catch, end

        % chi
        try
            [~, ~, mean_autocorr_seg, lm_ac] = autocorrelationDecay(ST_seg, lenST_seg, 1/3);
            seg_chi(seg) = -lm_ac.Coefficients.Estimate(2);
        catch, end

    catch
        % fit failed — leave NaN
    end
end

% --- 5. Count avalanches per bin ---
n_total    = zeros(1, n_bins);
n_powerlaw = zeros(1, n_bins);
for b = 1:n_bins
    t1_bin = bin_edges(b);
    t2_bin = bin_edges(b+1);
    in_bin = aval_starts >= t1_bin & aval_stops <= t2_bin;
    n_total(b) = sum(in_bin);

    bin_center = bin_centers(b);
    seg_idx = find(all_segments(:,1) <= bin_center & all_segments(:,2) >= bin_center, 1);
    if isempty(seg_idx) || isnan(seg_xminS(seg_idx)) || isnan(seg_xmaxS(seg_idx)) || isnan(seg_xminT(seg_idx)) || isnan(seg_xmaxT(seg_idx))
        continue
    end
    xminS = seg_xminS(seg_idx);
    xmaxS = seg_xmaxS(seg_idx);
    xminT = seg_xminT(seg_idx);
    xmaxT = seg_xmaxT(seg_idx);
    if xminS > xmaxS
        continue
    end
    if xminT > xmaxT
        continue
    end
    S_bin = S_all(in_bin);
    T_bin = T_all(in_bin);
    is_powerlaw = S_bin >= xminS & S_bin <= xmaxS & T_bin >= xminT & T_bin <= xmaxT;
    n_powerlaw(b) = sum(is_powerlaw);
end

% --- 6. Plot ---
fig_title = sprintf('Criticality over time — %s — %s', region, state);
[fig, axs] = makeFigure(fig_title, fig_title, [2,1]);
ax     = axs(1);   % main plot (top)
ax_stats = axs(2); % statistics (bottom)

% colored background bands per state (only when state='all')
if state == "all"
    state_segments = [];
    for s = 1:length(opt.states)
        try
            [~,~,s_index,~] = R.arrayInd(opt.states(s), region);
            st = R.state.times{s_index};
            state_segments = [state_segments; st, repmat(s, size(st,1), 1)];
        catch
            warning('plotStateCriticalityOverTime:missingState', 'State "%s" not found, skipping.', opt.states(s));
        end
    end
    state_segments = sortrows(state_segments, 1);
    for s = 1:length(opt.states)
        color = opt.colors(s,:);
        segs_s = state_segments(state_segments(:,3)==s, 1:2);
        for k = 1:size(segs_s,1)
            if k == 1, vis = 'on'; else, vis = 'off'; end
            patch(ax, [segs_s(k,1) segs_s(k,2) segs_s(k,2) segs_s(k,1)], ...
                [0 0 max(n_total)*1.2 max(n_total)*1.2], color, ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                'DisplayName', char(opt.states(s)), ...
                'HandleVisibility', vis);
        end
    end
end

plot(ax, bin_centers, n_total,    'k-',  'LineWidth', 2, 'DisplayName', 'Total avalanches');
plot(ax, bin_centers, n_powerlaw, 'r--', 'LineWidth', 2, 'DisplayName', 'Power-law avalanches');

xlabel(ax, 'Time (s)');
ylabel(ax, 'Number of avalanches per bin');
xlim(ax, [t_start, t_end]);
if max(n_total) > 0
    ylim(ax, [0, max(n_total)*1.2]);
end
legend(ax, 'show');

% --- 7. Plot statistics per segment ---
seg_centers = (all_segments(:,1) + all_segments(:,2)) / 2;
valid = ~isnan(seg_alpha);

stats = {seg_alpha, seg_beta, seg_gamma_exp, seg_gamma_area, seg_gamma_shape, seg_chi};
labels = {'\alpha', '\beta', '\gamma_{exp}', '\gamma_{area}', '\gamma_{shape}', '\chi'};
colors_stats = lines(length(stats));

%display(seg_alpha(1:4));
%display(valid);

hold(ax_stats, 'on');
for k = 1:length(stats)
    v = stats{k};
    plot(ax_stats, seg_centers(valid), v(valid), 'o-', ...
        'Color', colors_stats(k,:), 'LineWidth', 1.5, ...
        'DisplayName', labels{k});
end
xlabel(ax_stats, 'Time (s)');
ylabel(ax_stats, 'Exponent value');
xlim(ax_stats, [t_start, t_end]);
legend(ax_stats, 'show');

hold(ax_stats, 'off');
hold(ax, 'off');
end



%intervals_all = Restrict(intervals_all, state_times, shift=true);
%S_all = Restrict([intervals_all(:,1),S_all], state_times, shift=true);
%S_all = S_all(:,2);
%T_all = Restrict([intervals_all(:,1),T_all], state_times, shift=true);
%T_all = T_all(:,2);

%aval_starts = intervals_all(:,1);
%aval_stops  = intervals_all(:,2);