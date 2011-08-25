function res = results(fil,fir)
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
    res.freqs = (0:res.N/2-1)'*fil.Fs/res.N;
    res.lambda = 1/res.N*diag(res.e'*res.e);
    
    res.E = 1/sqrt(res.N)*fft(res.e-repmat(mean(res.e),res.N,1));
    res.E(res.N/2+1:end,:) = [];
    res.See = conj(res.E).*res.E;
end