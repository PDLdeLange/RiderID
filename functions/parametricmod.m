function [mod] = parametricmod(fir,bike,flag)
% Script for parametric fitting on the impulse response data.

    % Initial parameters based on optimal control theory:
    K = optimal(bike);
    % Preallocating some field names (prevents trouble)
    mod = struct('G',[],'K',[],'X',[],'X0',[],'z',[],'y',[],'C',[]); 
    % Initial conditions for the different model structurs
    mod(1).X0 = [  -33.1474   -1.8789  -22.1774    4.5624  -73.0965    0.7492 0];
    mod(2).X0 = [  -33.4989         0  -22.2414    4.5369  -74.1971    0.7441 0];
    mod(3).X0 = [  -42.4170         0  -14.5483         0  -61.3895    0.7441 0];
    mod(4).X0 = [  -42.4170         0  -14.5483   11.6824  -61.3895         0 0];
    mod(5).X0 = [  -32.7267   -1.3957  -28.6244    1.3743  -76.5304    2.7280 0.06];
    mod(6).X0 = [  -31.0394         0  -28.6244    1.3743  -76.5304    2.7280 0.06];
    mod(7).X0 = [  -28.8509         0  -29.9154         0  -79.1461    2.8067 0.06];
    
% RANDOM INITIAL PARAMETER SEARCH METHOD -> BRUTE FORCE APPROACH
%     for k = 1:500
%     mod(4).X0{k} = 50*(randn(1,6)-0.5);
%     en = errorfunc(ones(size(mod(3).X0{k})),mod(3).X0{k},1,fir,bike,3);
%     V(k) = en'*en;
%     end
%     keyboard;
    
    % For each model structure
    for i = 1:length(mod);
        % Bike
        mod(i).G.yu =  bike(3:4,2);
        mod(i).G.yw =  bike(3:4,3);
        mod(i).G.zu = -bike(3,2);
        mod(i).G.zw = -bike(3,3);
        % Normalizing variables using initial parameter conditions.
        X0 = mod(i).X0;
        X0n = ones(size(X0));
        theta0 = X0(logical(X0));  %#ok<*AGROW>
        theta0n = X0n(logical(X0));  %#ok<*AGROW>
        e0  = norm(errorfunc(theta0n,X0,1,fir,i,mod(i)));
        % Optimize model using the LSQNONLIN algorithm if flag is set to 1
        if flag == 1         
            % Optimize!
            [thetan,resnorm,en,exitflag,output,~,Jn] = ...
                lsqnonlin(@(thetan)errorfunc(thetan,X0,e0,fir,i,mod(i)),theta0n);    
            % Unnormalizing output
            Xn = zeros(size(X0)); Xn(logical(X0)) = thetan;
            X = Xn.*X0;
            theta = thetan.*theta0;
            e = en.*e0;
            J = Jn./repmat(theta0,size(Jn,1),1).*e0;
            % Optimization output
            N = length(e);
            mod(i).covP = abs(1/N*(e'*e)*inv(J'*J)); %#ok<MINV>
            mod(i).covPn = full(mod(i).covP)./(theta'*theta);
            mod(i).sem = full(sqrt(diag(inv(J'*J))*sum(e.^2)/N));
            mod(i).resnorm = resnorm;
            mod(i).exitflag = exitflag;
            mod(i).output = output;
            mod(i).sel = logical(X0);
        else % dont optimize
            X = X0n.*X0;
            disp('Skipping optimiziation, turn flag to 1 to optimize');
        end
        % Model output
        mod(i).i = i;
        mod(i) = riderfunc(X,tf('s'),i,mod(i)); 
        
        u = zeros(size(fir.tau)); u(1) = 1/(fir.tau(2)-fir.tau(1));
        g_est = lsim(mod(i).y,u,fir.tau-fir.tau(1));
        g = [fir.g(fir.m+1:end,:);zeros(fir.m,2)];
        mod(i).vaf = vaf(g,g_est);
    end  
end

% Error definition
function en = errorfunc(thetan,X0,e0,fir,i,mod)
    Xn = zeros(size(X0)); Xn(logical(X0)) = thetan;
    X = Xn.*X0; % Renormalizing
    mod = riderfunc(X,tf('s'),i,mod); % X,s,i,mod
%     ymod = impulse(mod.y,fir.tau-fir.tau(1));
    u = zeros(size(fir.tau)); u(1) = 1/(fir.tau(2)-fir.tau(1));
    ymod = lsim(mod(1).y,u,fir.tau-fir.tau(1));
    gmod = ymod;
    g = [fir.g(fir.m+1:end,:);zeros(fir.m,2)];
    etmp = (g - gmod)./repmat(std(g),size(g,1),1);
    en = [etmp(:,1); flipud(etmp(:,2))]/e0;
    disp(sum(en.^2));
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


