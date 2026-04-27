function profile_average = AvalAverageSizeTimeDependent(profiles,lifetimes,opt)
% average profile for avalanches of the same lifetime
% profiles : cell array, every element is an avalanche profile (column vector)
% lifetimes : lifetime of every avalanche, optional (can be deduced from 'profiles')
% min_lifetime : option, discard avalanches with smaller lifetime

arguments
  profiles
  lifetimes = []
  opt.min_lifetime = 4
end

if isempty(lifetimes)
  lifetimes = cellfun(@(x) size(x,1), profiles);
end

profile_average = cell(max(lifetimes),1);
for i = opt.min_lifetime : numel(profile_average)
  these_profiles = [profiles{lifetimes == i}];
  profile_average{i} = mean(these_profiles,2);
end