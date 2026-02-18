function plotRawRegionAutocorrDecay(regions, T_ref, opt)
%wrapper rawRegionAutocorrDecay and plotAutocorrelationDecay
arguments
    regions
    T_ref
    opt.states = {'all'}
    opt.regions = [-1]
end
    [~, autocorr, mean_autocorr, lm] = rawRegionAutocorrDecay(regions, T_ref, states=opt.states, regions=opt.regions);
    plotAutocorrelationDecay(autocorr, mean_autocorr, lm, T_ref);
end