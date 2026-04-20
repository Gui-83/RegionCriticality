 function plotCriticalityOverTime(R, region, state, bin_size, window, opt)
% plotCriticalityOverTime Plot number of avalanches and power-law avalanches over time
% with colored background bands for each state
%
% arguments:
%   R           - regions object
%   region      - string, brain region
%   bin_size    - double, bin size in seconds
%
% name-value arguments:
%   interval    - (1,2) double, [start, stop] of session (default: full session)
%   states      - string array of states (default: ["sws","rem","other"])
%   colors      - (n_states,3) double, RGB colors per state
%   min_avals   - minimum avalanches in a state segment to attempt power law fit
%   p_threshold - p-value threshold for swipeMLEDiscretePowerLawBounds
%   n_swipe     - number of steps for swipe (passed to swipeMLEDiscretePowerLawBounds)
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
intervals_all = R.avalIntervals(state, region);
S_all = R.avalSizes(state, region)*window;
T_all = round((intervals_all(:,2) - intervals_all(:,1)) / window); % durations in bins

aval_starts = intervals_all(:,1);
aval_stops  = intervals_all(:,2);

% session interval
if any(isnan(opt.interval))
    t_start = min(aval_starts);
    t_end   = max(aval_stops);
else
    t_start = opt.interval(1);
    t_end   = opt.interval(2);
end

% --- 2. Build time bins ---
bin_edges  = t_start : bin_size : t_end;
if bin_edges(end) < t_end
    bin_edges(end+1) = t_end;
end
n_bins = length(bin_edges) - 1;
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

% --- 3. Compute [xmin, xmax] per state segment ---
% collect all state segments in chronological order
state_segments = []; % each row: [t1, t2, s_index]
for s = 1:length(opt.states)
        sstate = opt.states(s);
        try
            [~,~,s_index,~] = R.arrayInd(sstate, region);
            state_times = R.state.times{s_index};
            for k = 1:size(state_times,1)
                state_segments = [state_segments; state_times(k,1), state_times(k,2), s];
            end
        catch
            warning('plotCriticalityOverTime:missingState', 'State "%s" not found, skipping.', sstate);
        end
end

if isempty(opt.zones)
    all_segments = state_segments;
else
    all_segments = [opt.zones, zeros(size(opt.zones,1),1)]; % s_index = 0 (pas de state associé)
end 
% sort segments chronologically
all_segments = sortrows(all_segments, 1);


%Pour ne compacter qu'un state :
if state ~= "all"
    [~,~,s_index,~] = R.arrayInd(state, region);
    state_times = R.state.times{s_index};
    
    % shifter les starts et stops des avalanches
    [aval_starts_shifted, ind1] = Restrict(aval_starts, state_times, shift=true);
    [aval_stops_shifted,  ind2] = Restrict(aval_stops,  state_times, shift=true);
    ind = intersect(ind1, ind2);
    %aval_starts = aval_starts_shifted(ismember(ind1, ind));
    %aval_stops  = aval_stops_shifted(ismember(ind2, ind));
    aval_starts = aval_starts(ismember((1:length(ind1))', find(ismember(ind1,ind))));
    aval_stops  = aval_stops(ismember((1:length(ind2))', find(ismember(ind2,ind))));
    S_all = S_all(ind);
    T_all = T_all(ind);
    
    % shifter all_segments et state_segments de la même façon
    t_start = min(aval_starts);
    t_end = max(aval_stops);
    if isempty(opt.zones)
        all_segments = [t_start, t_end, 0];
        state_segments = [t_start, t_end, 0];
    else
        all_segments = [opt.zones, zeros(size(opt.zones,1),1)]; % s_index = 0 (pas de state associé)
        state_segments = [];
    end 
    
    % reconstruire les bins sur le temps shifté
    bin_edges   = t_start : bin_size : t_end;
    if bin_edges(end) < t_end
        bin_edges(end+1) = t_end;
    end
    n_bins = length(bin_edges) - 1;
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
end


% for each segment, fit power law on avalanches within it
n_seg = size(all_segments, 1);
seg_xmin = NaN(n_seg, 1);
seg_xmax = NaN(n_seg, 1);
for seg = 1:n_seg
    t1 = all_segments(seg, 1);
    t2 = all_segments(seg, 2);
    in_seg = aval_starts >= t1 & aval_stops <= t2;
    S_seg = S_all(in_seg);
    T_seg = T_all(in_seg);
    if sum(in_seg) < opt.min_avals
        continue
    end
    try
        [~, xmin_S, xmax_S, ~] = swipeMLEDiscretePowerLawBounds(S_seg, opt.p_threshold, opt.n_swipe, 0.1, 3, false);
        [~, xmin_T, xmax_T, ~] = swipeMLEDiscretePowerLawBounds(T_seg, opt.p_threshold, opt.n_swipe, 1,   3, false);
        seg_xmin(seg) = max(xmin_S, xmin_T);
        seg_xmax(seg) = min(xmax_S, xmax_T);
    catch
        % not enough data or fit failed — leave NaN
    end
end

% --- 4. For each bin, count total and power-law avalanches ---
n_total   = zeros(1, n_bins);
n_powerlaw = zeros(1, n_bins);

for b = 1:n_bins
    t1_bin = bin_edges(b);
    t2_bin = bin_edges(b+1);
    in_bin = aval_starts >= t1_bin & aval_stops <= t2_bin;
    n_total(b) = sum(in_bin);

    % find which segment this bin belongs to
    bin_center = bin_centers(b);
    seg_idx = find(all_segments(:,1) <= bin_center & all_segments(:,2) >= bin_center, 1);
    if isempty(seg_idx) || isnan(seg_xmin(seg_idx)) || isnan(seg_xmax(seg_idx))
        continue
    end
    xmin = seg_xmin(seg_idx);
    xmax = seg_xmax(seg_idx);
    if xmin > xmax
        continue
    end
    S_bin = S_all(in_bin);
    T_bin = T_all(in_bin);
    is_powerlaw = S_bin >= xmin & S_bin <= xmax & T_bin >= xmin & T_bin <= xmax;
    n_powerlaw(b) = sum(is_powerlaw);
end

% --- 5. Plot ---
fig_title = sprintf('Criticality over time — %s', region);
[~, ax] = makeFigure(fig_title, fig_title, [1,1]);

% colored background bands per state
if state == "all"
    for s = 1:length(opt.states)
        color = opt.colors(s,:);
        seg_mask = state_segments(:,3) == s;
        segs_s = state_segments(seg_mask, 1:2);
        for k = 1:size(segs_s,1)
                if k == 1
                    vis = 'on';
                else
                    vis = 'off';
                end
                patch(ax, [segs_s(k,1) segs_s(k,2) segs_s(k,2) segs_s(k,1)], ...
                    [0 0 max(n_total)*1.2 max(n_total)*1.2], color, ...
                    'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
                    'DisplayName', char(opt.states(s)), ...
                    'HandleVisibility', vis);
        end
    end
end

% curves
plot(ax, bin_centers, n_total,    'k-',  'LineWidth', 2, 'DisplayName', 'Total avalanches');
plot(ax, bin_centers, n_powerlaw, 'r--', 'LineWidth', 2, 'DisplayName', 'Power-law avalanches');

xlabel(ax, 'Time (s)');
ylabel(ax, 'Number of avalanches per bin');
xlim(ax, [t_start, t_end]);
ylim(ax, [0, max(n_total)*1.2]);
legend(ax, 'show');
hold(ax, 'off');
end