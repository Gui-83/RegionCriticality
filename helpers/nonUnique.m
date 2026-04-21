function y = nonUnique(x)
  % remove outliers that appear only once
  [~,~,idx] = unique(x(:));
  counts = accumarray(idx,1);
  y = x(counts(idx) > 1);
end