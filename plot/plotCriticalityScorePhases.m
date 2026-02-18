function plotCriticalityScorePhases(Ids, States, Phases, G)
    %the inputs are the output of rawRegionSwipePhases
    n = length(Ids);
    figure;
    row = 3;
    column = ceil(n/row);
    for i = 1:n
        subplot(row, column, i);
        ids = Ids{i};
        states = States{i};
        score = 1./var(G{i}, 0, 3); %this is the criticality score (i.e the closer the 3 estimation of gamma the better the criticality)
        plotCriticalityScore(ids, states, score);
        title(sprintf('Criticality score %s', Phases{i}));
    end
end

