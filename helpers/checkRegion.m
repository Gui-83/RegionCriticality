function flag = checkRegion(region, regions)
    %deprecated function = raise an error if -1 is given as an argument, it
    %previouly meant 'take all the regions' (yes I know 0 did that very well
    %but I realized that too late and kept -1 to be backward compatible)
    % change 0 by the number that represent 'all the regions'

    if any(regions == -1)
        error('regions parameter must be provided. See README.md - Compute overall analysis from an instance of regions.');
    end


    if ismember(-1, regions)
        flag = ~(region.id == 0);
    else
        flag = ismember(region.id, regions);
    end
end