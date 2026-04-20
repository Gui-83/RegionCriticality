function [p,h] = StatTests(stat)

%x = [];
%n_conditions=length(group);
%for c= 1:n_conditions
    %x = [x;stat.(group(c))(:)];
%end
%[p,h,values,errors] = ANOVATests(x,group,'paired', true, 'parametric', false);

x = [stat.ISR(:), stat.other(:), stat.rem(:), stat.sws(:), stat.sws_non_ISR(:)];
x = x(~any(isnan(x), 2), :);  % retirer les lignes avec au moins un NaN
groups = [];
[p, h] = ANOVATests(x, groups, 'parametric', false, 'paired', false);