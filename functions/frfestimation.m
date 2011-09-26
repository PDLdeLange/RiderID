function frf = frfestimation(dat,set)

    t = dat.t;
    f = linspace(-1,1,set.N)'*set.Fs/2;
    N = set.N; frf.N = N;

    window = hann(N);

    y = dat.y.*repmat(window,1,2); %y = y-repmat(mean(y),N,1);
    w = dat.w.*window; %w = w-repmat(mean(w),N,1);

    frf.y = y;
    frf.w = w;
    
    Y = fftshift(1/sqrt(N)*fft(y),1); % FFT(y(t))
    W = fftshift(1/sqrt(N)*fft(w),1); % FFT(w(t))    

    Syy = conj(Y).*Y;
    Syw = conj(Y).*repmat(W,1,2);
    Sww = conj(repmat(W,1,2)).*repmat(W,1,2);
    
    G = Syw./Sww;
    
    m = 2^3+1; 
    
    mSyy = freqwinavg(Syy,f,m);
    mSyw = freqwinavg(Syw,f,m);
    mSww = freqwinavg(Sww,f,m);
    
    mG = mSyw./mSww;
    
    mSvv = mSyy - abs(mG).^2.*mSww;

    % Frequency Results
    frf.f = f;
    frf.W = W;
    frf.Y = Y;
    frf.Syy = mSyy;
    frf.Syw = mSyw;
    frf.Sww = mSww;
    frf.Svv = mSvv;
    frf.G = mG;
    
    % Impulse Response Results
    tau = dat.t-dat.t(ceil(N/2));
    gtmp = real(ifft(sqrt(N)*ifftshift(mG,1)));
    Rvtmp = real(ifft(sqrt(N)*ifftshift(mSvv,1)));
    
    g = zeros(size(gtmp));
    g(ceil(N/2):end,:) = flipud(gtmp(ceil(N/2):end,:)); g(tau==0,:) = 0;
    Rv = zeros(size(Rvtmp));
    Rv(ceil(N/2):end,:) = flipud(Rvtmp(ceil(N/2):end,:));
    
    frf.tau = tau;
    frf.g = g;    
    
    lambda = Rv(tau==0,:);
    H2 = mSvv./repmat(lambda,N,1);
    frf.lambda = lambda;
    frf.H2 = H2;
    
    % Digital bandpass filter
    L = zeros(size(y));
    L(f>=set.f1 & f<=set.f2,:) = 1;
    L(f<=-set.f1 & f>=-set.f2,:) = 1;
    frf.L = L;
   
   
end
    

