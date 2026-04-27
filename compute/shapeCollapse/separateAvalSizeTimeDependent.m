function [profiles,lifetimes] = separateAvalSizeTimeDependent(size_t)
% separate avalanche profiles, useful for shape collapse
% size_t : vector of avalanche profiles separated by zeros

if isempty(size_t)
  profiles = {};
  return
end
if size_t(1) ~= 0
  size_t = [0;size_t];
end
if size_t(end) ~= 0
  size_t = [size_t;0];
end

zerosPos = find(size_t == 0);
n = length(zerosPos);
profiles = cell(n-1,1);
lifetimes = zeros(n-1,1);
for i = 2 : n
  profiles{i-1} = size_t(zerosPos(i-1)+1 : zerosPos(i)-1);
  lifetimes(i-1) = zerosPos(i) - 1 - zerosPos(i-1);
end