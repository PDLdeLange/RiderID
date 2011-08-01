function res = results(fil,fir,mod)
% Testing the fir filter by simulating it using the measured input f.

    % Simulation by convolution of the FIR model and input force f.]
    res = fil;
    for i = 1:size(fil.y,2);
        res.y(:,i) = conv(fil.f,fir.g(:,i),'same')/fil.Fs;
    end

    % Error is difference between model and simulation output
    res.N = size(res.y,1);
    res.e = fil.y - res.y; 
    res.legend = {'\phi','\delta'};
    res.ylim = 4/3*max(abs(res.y))'*([-1 1]+1/4);
    
    
    res.E = 1/sqrt(res.N)*fft(res.e-repmat(mean(res.e),res.N,1));
    res.See = conj(res.E).*res.E;
end