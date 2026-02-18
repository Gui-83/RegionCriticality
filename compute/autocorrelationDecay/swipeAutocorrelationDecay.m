function chi = swipeAutocorrelationDecay(ST, opt)
%swipe the reference time parameter T_ref to find the best fit of slow
%autocorelation decay. See Tian et .al Theoretical foundations of studying
%criticality in the brain

%The best chi estimation is the one with the largest r squared
arguments
    ST
    opt.step = 10
    opt.rSquaredMin = 0.9
    opt.interpn = 1000
end
    T = linspace(0, 1, opt.step);
    chi = zeros(1, opt.step);
    f = waitbar(0 / opt.step, sprintf('Sending jobs (%d/%d)', 0, opt.step));
    future(1:opt.step) = parallel.FevalFuture;
    for i = 1:opt.step
        T_ref = T(i);
        future(i) = parfeval(@autocorrTask, 1, ST, T_ref, opt.rSquaredMin, opt.interpn);
        waitbar(i / opt.step, f, sprintf('Sending jobs (%d/%d)', i, opt.step));
    end
    
    waitbar(0, f, 'Pending...');

    for i = 1:opt.step
        [completedIdx, result] = fetchNext(future);
        chi(completedIdx) = result;
        waitbar(i / opt.step, f, sprintf('Computing (%d/%d)', i, opt.step));
    end
    close(f);
end

function chi = autocorrTask(ST, T_ref, rSquaredMin, interpn)
    [~, ~, ~, lm] = autocorrelationDecay(ST, T_ref, rSquaredMin=rSquaredMin, interpn=interpn);
    chi = -lm.Coefficients.Estimate(2);
end