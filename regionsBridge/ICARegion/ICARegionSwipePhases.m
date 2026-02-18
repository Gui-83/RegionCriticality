function [Ids, States, Phases, G] = ICARegionSwipePhases(session, opt)
arguments
  session
  opt.significanceLevel = 0.2
  opt.dicoStep = 10
  opt.min_decade = 1
  opt.base = 1.5
  opt.threshold = 4
end
    phases = loadEvents(session);
    states = ["awake"; "rem"; "sws"];
    n = length(phases);
    G = cell(1, n);
    Ids = cell(1, n);
    States = cell(1, n);
    Phases = cell(1, n);
    f = waitbar(0,sprintf('Computing (%d/%d)', 0, n));
    for i = 1:n
        phase = phases(i);
        R = regions(session, phases=phase, states=states);
        R = R.loadSpikes();
        R = R.computeAvalanches();
        gammas = ICARegionSwipe(R, significanceLevel=opt.significanceLevel, dicoStep=opt.dicoStep, min_decade=opt.min_decade, base=opt.base, threshold=opt.threshold);
        G{i} = gammas;
        Ids{i} = R.ids;
        States{i} = R.states;
        Phases{i} = phase;
        waitbar(i/n, f, sprintf('Computing (%d/%d)', i, n));
    end
    close(f);
end