function timeDependendentSize = catICAAvalSizeTimeDependent(region)
    timeDependendentSize = [];
    for i = 1:numel(region.brain_array)
        brain = region.brain_array(i);
        timeDependendentSizeBrain = brain.aval_timeDependendentSize;
        if timeDependendentSizeBrain(1) ~= 0
            timeDependendentSizeBrain = [0;timeDependendentSizeBrain];
        end
        timeDependendentSize = [timeDependendentSize; timeDependendentSizeBrain];
    end
end