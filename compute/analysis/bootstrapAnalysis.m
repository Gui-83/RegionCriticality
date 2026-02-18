function [S, T, ST, Ealpha, Ebeta, Egamma1, Egamma2, Egamma3] = bootstrapAnalysis(S, T, ST, opt)
%use bootstrap to compute the distribution of the different gamma estimator
arguments
  S %size
  T %lifetime
  ST %time dependent size (see separateAvalSizeTimeDependent function)
  opt.significanceLevel = 0.05 %significant level for power law fit
  opt.dicoStep = 10 %dicotomie step for power law fit
  opt.min_decade = 1 %minimum decade for power law fit range
  opt.base = 1.5 %see powerLawFit.md
  opt.threshold = 4 %minimum lifetime of avalanche to consider in shape collapse estimation
  opt.samples = 0.05 %sample percentage of the raw data fot each sub sampling
  opt.bootstrap = 1000 %number of subsample to compute
end
    samples = floor(length(S)*opt.samples);
    [BS, BT, BST] = bootstrapping(S, T, ST, samples, opt.bootstrap); %sample op.bootstrap times the data

    %those are the list that will contain the result of the bootstrap
    Ealpha = zeros(1, opt.bootstrap);
    Ebeta = zeros(1, opt.bootstrap);
    Egamma1 = zeros(1, opt.bootstrap);
    Egamma2 = zeros(1, opt.bootstrap);
    Egamma3 = zeros(1, opt.bootstrap);
    
    f = waitbar(0 / opt.bootstrap, sprintf('Sending jobs (%d/%d)', 0, opt.bootstrap));

    future(1:opt.bootstrap) = parallel.FevalFuture;
    for i = 1:opt.bootstrap
        S = BS(:,i);
        T = BT(:,i);
        ST = BST{i};
        future(i) = parfeval(@boostraptask, 5, S, T, ST, opt);
        waitbar(i / opt.bootstrap, f, sprintf('Sending jobs (%d/%d)', i, opt.bootstrap));
    end
    
    waitbar(0, f, 'Pending...');

    for i = 1:opt.bootstrap
        [completedIdx, alpha, beta, gamma1, gamma2, gamma3] = fetchNext(future);
        Ealpha(completedIdx) = alpha;
        Ebeta(completedIdx) = beta;
        Egamma1(completedIdx) = gamma1;
        Egamma2(completedIdx) = gamma2;
        Egamma3(completedIdx) = gamma3;
        waitbar(i / opt.bootstrap, f, sprintf('Computing (%d/%d)', i, opt.bootstrap));
    end

    close(f)


    indices = ~isnan(Ealpha); % remove unsuccessful fit
    Ealpha = Ealpha(indices);
    Ebeta = Ebeta(indices);
    Egamma1 = Egamma1(indices);
    Egamma2 = Egamma2(indices);
    Egamma3 = Egamma3(indices);
    S = BS(:,indices);
    T = BT(:,indices);
    ST = BST(indices);
end

function [alpha, beta, area_gam, gam, shape_gam] = boostraptask(S, T, ST, opt)
    ST_lenght = cellfun(@length, ST);
    A = getArea(T, S);
    %we only keep the relevent part of the analysis function, meaning the apart that will let use compute the three estimator of gamma
    [alpha, ~, ~, ~, beta, ~, ~, ~, Alm, gam, ~, ~, ~, shape_gam] = analysis(S, T, A, ST, ST_lenght,ones(1, 3), significanceLevel=opt.significanceLevel, dicoStep=opt.dicoStep, min_decade=opt.min_decade, base=opt.base, threshold=opt.threshold, verbose=false);
    if alpha > 0 && beta > 0
        area_gam = Alm.Coefficients.Estimate(2);
    else
        [alpha, beta, area_gam, gam, shape_gam] = deal(NaN);
    end
end

