function plotCriticalityScore(ids, states, score, opt)
arguments
    ids
    states
    score
    opt.logscale = false;
end
    imagesc(score);
    yticks(1:length(ids));
    xticks(1:length(states));
    yticklabels(ids);
    xticklabels(states)
    colorbar;
    if opt.logscale
        set(gca,'ColorScale','log');
    end
end

