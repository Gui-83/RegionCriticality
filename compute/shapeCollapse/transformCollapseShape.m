function [x,S,T] = transformCollapseShape(shape,m,min_lifetime)
% interpolate shapes so that they all have length m
% x : linear time between 0 and 1

arguments
  shape
  m = 100
  min_lifetime = 1
end

shape = shape(cellfun(@numel, shape) >= min_lifetime);

n = numel(shape);
x = linspace(0,1,m);
S = zeros(m,n);
T = zeros(n,1);

for i = 1 : n
  s = shape{i};
  t = length(s);
  z = linspace(0,t,t);
  S(:,i) = interp1(z,s,x*t);
  T(i) = t;
end