function fil = davisfilter(raw,wn)
% Filters the data from the Davis Experiments. It removes some of the high
% frequency content on the state signals (y) and selects a range of data.

    % Copy data
    fil = raw;

    % Remove High frequency content
    [b,a] = butter(8,wn);
    fil.y = filtfilt(b,a,fil.y);

end

