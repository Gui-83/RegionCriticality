function [uniqueST_lengths, autocorr, mean_autocorr, lm] = autocorrelationDecay(ST, lengthST, T_ref, opt)
arguments
    ST
    lengthST
    T_ref
    opt.rSquaredMin = 0.8
    opt.interpn = 200
end
    uniqueST_lengths = unique(lengthST);
    uniqueST_lengths = uniqueST_lengths(uniqueST_lengths > 10); 
    uniqueST_lengths = sort(uniqueST_lengths);
    m = length(uniqueST_lengths);
    autocorr = cell(1, m);
    
    z = linspace(0, 1, opt.interpn);
    sum_autocorr = zeros(1, opt.interpn);
    numValidAutocorr = 0;
    for i = 1:m
        len = uniqueST_lengths(i);
        sameLenSTindices = lengthST == len;
        if sum(sameLenSTindices) < 20
            continue
        end
        sameLenST = ST(sameLenSTindices);
        sameLenSTMatrix = cell2mat(sameLenST);
        sameLenAutocorr = zeros(1, len);
        
        k = max(round(len*T_ref)+1, 1);
        ref = sameLenSTMatrix(k,:);
        n = length(ref);
        mean_ref = mean(ref);
        cov_ref = sum((ref - mean_ref) .^ 2)/(n-1);
        for j = 1:len
            cursor = sameLenSTMatrix(j,:);
            mean_cursor = mean(cursor);
            sameLenAutocorr(j) = sum((ref - mean_ref) .* (cursor - mean_cursor))/(n-1)/cov_ref;
        end
        autocorr{i} = sameLenAutocorr;
        sum_autocorr = sum_autocorr + interp1(linspace(0, 1, len), sameLenAutocorr, z);
        numValidAutocorr = numValidAutocorr + 1;
    end
    mean_autocorr = sum_autocorr/numValidAutocorr;
    
    indices = (z>T_ref);
    mean_autocorr(mean_autocorr <= 0) = min(abs(mean_autocorr));
    log_mean = log(mean_autocorr);
    
    x = log(z(indices) - T_ref);
    y = log_mean(indices);
    for i = 0:floor(length(x)/2)
        headx = x(1:end-i);
        heady = y(1:end-i);
        p = polyfit(headx, heady, 1);

        y_pred = polyval(p, headx);
        SS_tot = sum((heady - mean(heady)).^2);
        SS_res = sum((heady - y_pred).^2);
        Rsquared = 1 - (SS_res / SS_tot);
        n = length(headx);
        k = 1;
        Rsquared_adj = 1 - (1 - Rsquared) * (n - 1) / (n - k - 1);

        lm.Coefficients.Estimate = flip(p);
        lm.Rsquared.Adjusted = Rsquared_adj;
        if lm.Rsquared.Adjusted > opt.rSquaredMin
            break
        end
    end
end

