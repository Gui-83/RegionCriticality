function y = nonUnique(x,n_min)
% remove outliers that appear only once

arguments
  x
  n_min = 1
end

[~,~,idx] = unique(x(:));
counts = accumarray(idx,1);
y = x(counts(idx) > n_min);