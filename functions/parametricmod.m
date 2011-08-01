function [sys,opt] = parametricmod(fir,bike,flag)
% Script for determining an optimal parametric fit of the impulse response.

    % Initial parameters based on optimal control theory:
    K = optimal(bike);
    
    X0 = [  -15.2102  -57.3110  -19.4106    0.6815  -29.9533    6.2050];
%     X0 = [   -6.3163  -47.7830  -13.4008   -0.8079  -22.5410    6.0228];
    
    X0n = ones(size(X0)); 
    e0  = norm(errorfunc(X0n,X0,1,fir,bike));
    
    % Optimize the model using the LSQNONLIN algorithm if flag is set to 1
    if flag == 1
        [Xn,resnorm,en,exitflag,output,~,Jn] = lsqnonlin(@(Xn)errorfunc(Xn,X0,e0,fir,bike),X0n);
    else 
        Xn = X0n;
        disp('Skipping optimiziation, turn flag to 1 to enable optimization');
    end
    
    % Unnormalizing output
    e =  en.*e0;
    X  = Xn.*X0;
    J  = Jn./repmat(X0,size(Jn,1),1).*e0;
    sem = full(sqrt(diag(inv(J'*J))*sum(e.^2)/length(e)));

    
    % Model output
    sys = riderfunc(X,bike);
    sys.X = X;
    sys.lqr.K = K;
    
    % Optimization output
    opt.X0 = X0;
    opt.Jn = Jn;
    opt.Xn = Xn;
    opt.resnorm = resnorm;
    opt.en = en;
    opt.exitflag = exitflag;
    opt.output = output;
    opt.sem = sem;

end

% Error definition
function en = errorfunc(Xn,X0,e0,fir,bike)
    X = Xn.*X0; % Renormalizing
    sys = riderfunc(X,bike);
    ymod = impulse(sys.y,fir.tau-fir.tau(1));
    gmod = ymod;
    g = [fir.g(fir.m+1:end,:);zeros(fir.m,2)];
    en = reshape(g - gmod,numel(g),1)/e0;
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


