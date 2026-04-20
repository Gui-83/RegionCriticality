function stats_list = BatchCriticalViolinExponents(sessions, region)
arguments
    sessions
    region (1,1) string = "nr"
end
states_list = ["all", "other", "rem", "sws", "sws"];
event = "slownr";
event2 = "exceptslownr";
event_list = [event, "all", "all", "all", event2];
%region = "nr";
window = 0.05;

stats_names = {'\alpha', '\beta', '\gamma_{exp}', '\gamma_{area}', '\gamma_{shape}', '\chi', '\gamma_{area} - \gamma_{exp}'};
n_stats = 7;
n_conditions = 5; % sws, other, rem, slownr

[lessessions,extra_args] = readBatchFile(sessions);
n_sessions = length(lessessions);
args = repmat({region, states_list, event_list, window}, n_sessions, 1);

[alpha, beta, expected_gamma, gam_area, gam_shape, chi, diffgam] = runBatch(sessions, @CriticalExponents, args);
% plot
alpha = reverseCellStruct(alpha);
beta = reverseCellStruct(beta);
expected_gamma = reverseCellStruct(expected_gamma);
gam_area = reverseCellStruct(gam_area);
gam_shape = reverseCellStruct(gam_shape);
chi = reverseCellStruct(chi);
diffgam = reverseCellStruct(diffgam);

%data = [results(1,:), results(2,:), results(3,:), results(4,:)];

stats_list = {alpha, beta, expected_gamma, gam_area, gam_shape, chi, diffgam};