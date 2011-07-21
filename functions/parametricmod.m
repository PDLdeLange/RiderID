function sys = parametricmod(fir,bike,flag)
% Script for determining an optimal parametric fit of the impulse response.

    % Initial parameters based on optimal control theory:
    K = optimal(bike);
    
    X0 = [  -15.2102  -57.3110  -19.4106    0.6815  -29.9533    6.2050];
%     X0 = [   -6.3163  -47.7830  -13.4008   -0.8079  -22.5410    6.0228];
    
    X0n = ones(size(X0)); 
    e0n = norm(errorfunc(X0n,X0,1,fir,bike));
    
    % Optimize the model using the LSQNONLIN algorithm if flag is set to 1
    if flag == 1
        Xn = lsqnonlin(@(Xn)errorfunc(Xn,X0,e0n,fir,bike),X0n);
    else 
        Xn = X0n;
        disp('Skipping optimiziation, turn flag to 1 to enable optimization');
    end
    
    % Evaluate model
    X  = Xn.*X0;    
    sys = riderfunc(X,bike);
    sys.X = X;
    sys.X0 = X0;
    sys.lqr.K = K;

end

% Error definition
function en = errorfunc(Xn,X0,e0n,fir,bike)
    Fs = 1/(fir.tau(2)-fir.tau(1));
    X = Xn.*X0; % Renormalizing
    sys = riderfunc(X,bike);
    ymod = impulse(sys.y(1,:),fir.tau-fir.tau(1));
    gmod = ymod/Fs;
    g = [fir.g(fir.m+1:end,1);zeros(fir.m,1)];
    en = (g - gmod)/e0n;
    disp(1/2*sum(en.^2));
end

% Optimal controller
function K = optimal(bike)

    % Set performance weightning Q and control effort R
    q1max = 0.1; Q33 = 1/q1max^2;
    T2max = 2.0; R11 = 1/T2max^2;  

    % Weighting matrices
    Q = zeros(4); Q(3,3) = Q33;
    R = zeros(1); R(1,1) = R11;

    % Controller
    [K,~,e] = lqr(bike(:,2),Q,R); disp(e);

end


