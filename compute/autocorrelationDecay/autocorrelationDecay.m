function [unique_lifetime, autocorr, mean_autocorr, lm] = autocorrelationDecay(profile, lifetime, T_ref, opt)
% profile : cell array, each element is an avalanche profile (column vector)
% lifetime : int, each element is the lifetime of an avalanche in 'profile'
% min_lifetime = 10 : discard avalanches with shorter lifetime
% min_samples = 20 : ignore lifetimes which have less than 'min_samples' avalanches

arguments
  profile
  lifetime
  T_ref
  opt.min_lifetime = 10
  opt.min_samples = 20
  opt.rSquaredMin = 0.8
  opt.interpn = 200
end

% apply 'min_lifetime' and 'min_samples' filtering
[n_samples,unique_lifetime] = groupcounts(lifetime);
valid = (n_samples >= opt.min_samples & unique_lifetime >= opt.min_lifetime);
unique_lifetime = unique_lifetime(valid);

m = length(unique_lifetime);
autocorr = cell(1,m);
    
z = linspace(0,1,opt.interpn);
sum_autocorr = zeros(1,opt.interpn);
for i = 1 : m
  len = unique_lifetime(i);
  these_profiles = [profile{lifetime == len}];
  this_autocorr = zeros(1,len);

  % extract metrics across avalanches at reference time
  t_ref = max(round(len*T_ref)+1,1);
  ref = these_profiles(t_ref,:);
  mean_ref = mean(ref);
  n = length(ref);
  cov_ref = sum((ref - mean_ref) .^ 2)/(n-1);

  % metrics at all times w.r.t. reference time
  for j = 1:len
    cursor = these_profiles(j,:);
    mean_cursor = mean(cursor);
    this_autocorr(j) = sum((ref - mean_ref) .* (cursor - mean_cursor)) / (n-1) / cov_ref;
  end

  autocorr{i} = this_autocorr;
  sum_autocorr = sum_autocorr + interp1(linspace(0,1,len),this_autocorr,z);
end
mean_autocorr = sum_autocorr / m;
    
indices = z > T_ref;
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

