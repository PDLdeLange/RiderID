 function [mod] = parametricmod(npm,dat,set)
% Script for parametric fitting on the impulse response data.

    % Initial parameters based on optimal control theory:
%     K = optimal(bike);
    bike = davisbike(set.v); % Bicycle model from Davis

  
    % Preallocating some field names (prevents trouble)
    mod = struct('G',[],'K',[],'X',[],'X0',[],'z',[],'y',[],'C',[]); 

    % Initial Parameters
    mod(01).X0 = [  -77.7975   -2.4529  -57.1454    1.9939    7.9663 -174.4105    6.9691    0.1959];      
    mod(02).X0 = [  -78.3373         0  -57.2744    1.9845    7.9175 -176.0056    6.9606    0.1975];      
    mod(03).X0 = [  -54.4049         0  -50.2570    0.8928         0 -147.7373    5.3292    0.1247];      
    mod(04).X0 = [  -51.8316         0  -43.3446         0         0 -131.1479    3.9576    0.1402];      
    mod(05).X0 = [  -36.7135         0  -32.9633         0         0  -89.2661    3.2486         0];      

    mod(06).X0 = [  -65.4008    3.0592  -48.0903   -0.0751    2.7997 -146.8299    5.0620    0.2618];      
    mod(07).X0 = [  -66.4620    2.2533  -48.5222         0    3.2212 -148.0136    5.1718    0.2647];      
    mod(08).X0 = [  -65.8637         0  -48.3804         0    3.2145 -146.4612    5.1698    0.2631];      
    mod(09).X0 = [  -57.6751         0  -49.5773         0         0 -145.6462    5.2120    0.2397];      

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
        e0  = norm(errorfunc(theta0n,X0,1,npm,i,mod(i),dat));
        % Optimize model using the LSQNONLIN algorithm if flag is set to 1
            % Optimize!
            [thetan,resnorm,en,exitflag,output,~,Jn] = ...
                lsqnonlin(@(thetan)errorfunc(thetan,X0,e0,npm,i,mod(i),dat),...
                  theta0n);
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

        % Model output
        mod(i).i = i;
        mod(i) = riderfunc(X,tf('s'),i,mod(i));
        
        delta_mod = lsim(mod(i).y(2),dat.w,dat.t);
        delta = npm.y(:,2);

        mod(i).vaf = vaf(delta,delta_mod);
        
    end  
end


% Optimal controller
function K = optimal(bike)
    % Set performance weightning Q and control effort R
    q1max = 0.1; Q33 = 1/q1max^2;
    T2max = 2.0; R11 = 1/T2max^2;  
    % Weighting matrices
%     Q = zeros(6); Q(3,3) = Q33; Q(5,5) = Q33;
    Q = 1*eye(6); Q(3,3) = Q33;
    R = zeros(1); R(1,1) = R11;
    % Controller
    [K,~,e] = lqr(bike(:,2),Q,R); disp(e);
end

% Simulation by convolution of the FIR model and input force f.]



% Error definition
function en = errorfunc(thetan,X0,e0,npm,i,mod,dat)
    Xn = zeros(size(X0)); Xn(logical(X0)) = thetan;
    X = Xn.*X0; % Renormalizing
    mod = riderfunc(X,tf('s'),i,mod); % X,s,i,mod
    delta_mod = lsim(mod.y(2),dat.w,dat.t);
    delta = npm.y(:,2);
    e = 1/npm.N*(delta - delta_mod);
    en = e/e0;
    disp(sum(en.^2));
end



% % Error definition
% function en = errorfunc(thetan,X0,e0,npm,i,mod,dat)
%     Xn = zeros(size(X0)); Xn(logical(X0)) = thetan;
%     X = Xn.*X0; % Renormalizing
%     mod = riderfunc(X,tf('s'),i,mod); % X,s,i,mod
%     u = zeros(size(npm.t(1:npm.m))); u(1) = 1/(npm.t(2)-npm.t(1));
%     ghat = lsim(mod.y,u,npm.t(1:npm.m));
%     g = npm.g(1:npm.m,:);
%     e = 1/npm.N*(g(:,2) - ghat(:,2));
%     en = e/e0;
%     disp(sum(en.^2));
% end
