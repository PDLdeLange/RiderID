function res = results(dat,fir,set)


% Testing the fir datter by simulating it using the measured input f.

    % Simulation by convolution of the FIR model and input force f.]
    for i = 1:size(dat.y,2);
        res.y(:,i) = conv(dat.w,[zeros(fir.m2+fir.m1,1); fir.g(:,i)],'same')/set.Fs;
    end
    
    % Error is difference between model and simulation output
    res.N = size(res.y,1);
    res.v = dat.y - res.y;
    res.legend = {'\phi','\delta'};
    res.ylim = 4/3*max(abs(res.y))'*([-1 1]+1/4);
    
    
    % Create prefilter L(omega)
    res.f1 =  0.1;
    res.f2 = 10.0;
    res.L = zeros(size(res.y));
    res.L(res.f>=res.f1 & res.f<=res.f2,:) = 1;
    res.L(res.f<=-res.f1 & res.f>=-res.f2,:) = 1;
    
    
    res.f = res.Fs/2*linspace(-1,1,res.N)';    
    res.V = fftshift(1/sqrt(res.N)*fft(res.v),1);
    res.Svv = conj(res.V).*res.V;
    
    rem = noisemod(res);    
    
    res.V_mod = squeeze(freqresp([rem.H(1) rem.H(2)],2*pi*res.f)).';
    res.v_mod = real(ifft(sqrt(res.N)*ifftshift(res.V_mod,1)));


    % Find noise model
    res.lambda = res.v_mod(1,:);
    res.h_mod = res.v_mod./repmat(res.lambda,res.N,1);
    res.H_mod = res.V_mod./repmat(res.lambda,res.N,1);

    

close all
loglog(res.f,abs(res.H_mod.*repmat(res.lambda,res.N,1)),':'); hold on;
loglog(res.f,abs(res.V))

close all
for i = 1:2
    figure(1);
    subplot(2,1,i);
    plot(res.v_mod(:,i),':'); hold on;
    plot(res.v(:,i));
    
    figure(2);
    subplot(2,1,i);
    hist([res.v_mod(:,i) res.v(:,i)],21)
end
    
    
end


function rem = noisemod(res)


    f = res.f; % frequency (Hz)
    V = res.V; % noise spectrum  v(j\omega)

    gamma = ones(size(V));

    X0 = [-35.5001 -48.2208 11.3043 38.187 -136.5499 -174.7491 6.6210 81.6086];

    X0n = ones(size(X0));
    e0 = error_freq(X0n,V,@V_modfunc,gamma,f,X0);

    lsqfunc = @(Xn)error_freq(Xn,V,@V_modfunc,gamma,f,X0);
    Xn = lsqnonlin(lsqfunc,X0n);
    X = Xn.*X0;
    % % function e = error_freq(X,H,H_modfunc,gamma,f)
    % lsqfunc = @(thetan)error_freq(thetan,mHarm_wb,@Harm_model,mCwx_wb,mf_wb,theta0,parms);
    % [thetan,resnorm,e,exitflag,output,ignore,Jn] = lsqnonlin(lsqfunc,theta0n);

    % He0 = He_modfunc(X0,2*pi*f);
    rem.f = f;
    rem.V = V;
    rem.V_mod = V_modfunc(X,f);
    rem.V0_mod = V_modfunc(X0,f);
    rem.X = X;

    s = tf('s');

    rem.H = [(X(1)+X(2)*s)/(s^2+X(3)*s+X(4))^2;
             (X(5)+X(6)*s)/(s^2+X(7)*s+X(8))^2;];

    rem.e0 = e0;
    rem.e = error_freq(Xn,V,@V_modfunc,gamma,f,X0);

end

% Error definition
function e = error_freq(Xn,V,V_modfunc,gamma,f,X0)
    X = Xn.*X0;
    V_mod = V_modfunc(X,f);
    e = reshape(sqrt(1./repmat(f,1,2)).*gamma.*...
        abs(log(V./V_mod)),numel(V),1);
end

function V_mod = V_modfunc(X,f)

    s = 1j*2*pi*f;
      
     V_mod = [(X(1)+X(2).*s)./(s.^2+X(3).*s+X(4)).^2 ... 
              (X(5)+X(6).*s)./(s.^2+X(7).*s+X(8)).^2;];
    

end
