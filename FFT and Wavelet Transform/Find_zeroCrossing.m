function zeroCrossings = Find_zeroCrossing(signal, min_gap)
    % Initialize the zero_crossings vector
    zeroCrossings = [];
    
    % Loop through the signal and find zero crossings with a minimum gap
    i = 1;
    while i < length(signal) - 1
        if (signal(i) >= 0 && signal(i+1) < 0) || (signal(i) < 0 && signal(i+1) >= 0)
            zeroCrossings = [zeroCrossings; (i + (i+1))/2];
            i = i + min_gap;
        else
            i = i + 1;
        end
    end
end