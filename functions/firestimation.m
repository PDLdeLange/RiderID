function fir = firestimation(fil,m,wn)

    % FIR finite impulse response function g(tau)
    N = size(fil.y,2);
    M = 2*m+1;

    % Timeshift vector tau
    tau = 1/fil.Fs*(-m:m)';

    % FIR windowing
    win = zeros(size(tau)); win(m+2:end) = 1;
    win(m+1+3/4*m+1:end) = 1/2*(1+cos((1:1/4*m)*4*pi/m));

    % Filter parameters
    [b, a] = butter(8,wn);

    % For every output y
    g = zeros(M,N); g_raw = zeros(M,N);
    for i = 1:N
        g_raw(:,i) = markovparm(fil.y(:,i),fil.f,m);
        g(:,i) = filtfilt(b,a,g_raw(:,i));
        g(:,i) = g(:,i).*win;
    end

    % Output
    fir.m = m;
    fir.tau = tau;
    fir.win = win;
    fir.g_raw = g_raw;
    fir.g = g;
    fir.legend = {'g_\phi','g_\delta'};

    
    % ARX model generation
    A(1:2,1:2,1) = eye(2);
    B = zeros(2,1,fir.m);
    B(1,1,:) = fir.g(fir.m+2:end,1); B(2,1,:) = fir.g(fir.m+2:end,2);
    fir.arx = idarx(A,B,1/fil.Fs);
    % fir.frd = idfrd(fir.arx);
    
end