function [bin_centers, n_total, n_powerlaw,n_proportion, t_start, t_end, all_segments, seg_alpha, seg_beta, seg_gamma_exp, seg_gamma_area, seg_gamma_shape, seg_chi] = StateCriticalityOverTime_sliding_window(R, region, state, bin_size, window, opt)

arguments
    R (1,1) regions
    region (1,1) string
    state (1,1) string = "all"
    bin_size (1,1) {mustBeNumeric, mustBePositive} = 0.5
    window (1,1) = 0.05
    opt.interval (1,2) double = [NaN, NaN]
    opt.p_threshold (1,1) double = 0.05
    opt.n_swipe (1,1) double = 10
    opt.st (1,1) double = 100
    opt.wd (1,1) double = 500
    opt.slownr (1,1) double = 0
end

% --- 1. Get all avalanches ---
if opt.slownr == 1
    ISRIntervals = R.eventIntervals("slownr");
    intervals_all = R.avalIntervals('all', region, 'restriction',ISRIntervals);
    S_all = R.avalSizes('all', region, 'restriction',ISRIntervals) * window;
else
    intervals_all = R.avalIntervals('all', region);
    S_all = R.avalSizes('all', region) * window;
end

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
a = 0:floor((t_end - opt.wd)/opt.st) - 1;
all_segments_zones = [a'*opt.st, a'*opt.st + opt.wd];
all_segments = [all_segments_zones, zeros(size(all_segments_zones,1),1)];

all_segments = sortrows(all_segments, 1);

% --- 4. Fit power law per segment ---
if state == "all"
    ind=[];
end
[seg_alpha, seg_beta, seg_gamma_exp, seg_xminS, seg_xmaxS, seg_xminT, seg_xmaxT,seg_gamma_area, seg_gamma_shape, seg_chi] = FitPoweLawPerSegment(R, region, S_all, T_all, aval_starts, aval_stops, all_segments, state, ind, p_threshold = opt.p_threshold);

% --- 5. Count avalanches per bin ---
n_total    = zeros(1, n_bins);
n_powerlaw = zeros(1, n_bins);
n_proportion = zeros(1, n_bins);
seg_centers = (all_segments(:,1) + all_segments(:,2)) / 2;

for b = 1:n_bins
    t1_bin = bin_edges(b);
    t2_bin = bin_edges(b+1);
    in_bin = aval_starts >= t1_bin & aval_stops <= t2_bin;
    n_total(b) = sum(in_bin);

    bin_center = bin_centers(b);
    %seg_idx = find(all_segments(:,1) <= bin_center & all_segments(:,2) >= bin_center, 1);
    [~, seg_idx] = min(abs(seg_centers - bin_center));
    if isnan(seg_alpha(seg_idx)) || isnan(seg_beta(seg_idx))
        continue
    end
    xminS = seg_xminS(seg_idx);
    xmaxS = seg_xmaxS(seg_idx);
    xminT = seg_xminT(seg_idx);
    xmaxT = seg_xmaxT(seg_idx);
    if xminS > xmaxS || xminT > xmaxT
        continue
    end
    S_bin = S_all(in_bin);
    T_bin = T_all(in_bin);
    is_powerlaw = S_bin >= xminS & S_bin <= xmaxS & T_bin >= xminT & T_bin <= xmaxT;
    n_powerlaw(b) = sum(is_powerlaw);
    if n_total(b) ~= 0
        n_proportion(b) = n_powerlaw(b)/n_total(b);
    end
end