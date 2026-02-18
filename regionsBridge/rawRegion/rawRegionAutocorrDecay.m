function [uniqueST_lengths, autocorr, mean_autocorr, lm] = rawRegionAutocorrDecay(regions, T_ref, opt)
%wrapper for autocorrelationDecay, extract data from region and feed it to
%the function. Returns the estimation of the slow decay of autorelation
%parameter
arguments
    regions
    T_ref
    opt.states = {'all'}
    opt.regions = [-1]
end
    timeDependendentSize = getRawRegionSizeTimeDependent(regions, states=opt.states, regions=opt.regions);
    [ST, lengthST] = separateAvalSizeTimeDependent(timeDependendentSize);
    [uniqueST_lengths, autocorr, mean_autocorr, lm] = autocorrelationDecay(ST, lengthST, T_ref);
end

