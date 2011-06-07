function [closestIndex] = lookupCDF(CDF, value)

    % Find and interpolate 'value' in 'CDF' and 'angles'
    minIndex = 2;
    maxIndex = length(CDF);
    midIndex = minIndex + ceil((maxIndex - minIndex)/2);

    while value ~= CDF(midIndex) && maxIndex >= minIndex

        midIndex = minIndex + ceil((maxIndex - minIndex)/2);
        if value > CDF(midIndex)
            minIndex = midIndex + 1;
        elseif value < CDF(midIndex)
            maxIndex = midIndex - 1;
        end

    end
    midIndex = minIndex + ceil((maxIndex - minIndex)/2);
    closestIndex = midIndex;
    
end