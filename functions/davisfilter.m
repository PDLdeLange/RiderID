function fil = davisfilter(raw,wn,range)
% Filters the data from the Davis Experiments. It removes some of the high
% frequency content on the state signals (y) and selects a range of data.

    % Copy data
    fil = raw;

    % Select a suitable range of data;
    fil.y = fil.y(range(1):range(2),:);
    fil.f = fil.f(range(1):range(2),:);
    fil.t = fil.t(range(1):range(2),:);

    % Remove High frequency content
    [b,a] = butter(8,wn);
    fil.y = filtfilt(b,a,fil.y);

end

