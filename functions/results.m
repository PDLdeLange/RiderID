function res = results(fil,fir,mod)
% Testing the fir filter by simulating it using the measured input f.

    % Simulation by convolution of the FIR model and input force f.]
    res = fil;
    for i = 1:size(fil.y,2);
        res.y(:,i) = conv(fil.f,fir.g(:,i),'same');
    end

    % Error is difference between model and simulation output
    res.e = fil.y - res.y; 
    res.legend = {'\phi','\delta'};
    res.ylim = 4/3*max(abs(res.y))'*([-1 1]+1/4);
        
end