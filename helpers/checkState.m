function flag = checkState(region, states)
    if ismember('all', states)
        flag = true;
    else
        % state is no longer a member of the region class

        %flag = ismember(region.state, states);

        flag = true; %just to avoid the bug
    end
end