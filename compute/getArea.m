function area = getArea(T,S)
% compute the Area given the lifetime and size

if isempty(T)
  area = [];
  return
end

area = accumarray(T,S,[],@mean);