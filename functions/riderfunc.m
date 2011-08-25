function [mod] = riderfunc(X,s,i,mod)

    % Declarations/limitations
    Gnm = 900/(s^2 + 2*.707*30*s + 900);
    delay = exp(-X(7)*s);

    mod.X = X(1:7);
    pid = [1;1/s;s];
    mod.C =  mod.X(1:6)*[pid zeros(3,1); zeros(3,1) pid];
    % Casewise multiple model structures
    switch i
        case 1; % PID controller
            mod.K = -mod.C*Gnm;
        case 2; % reduced PID type controller
            mod.K = -mod.C*Gnm;
        case 3; % reduced PID controller
            mod.K = -mod.C*Gnm;
        case 4; % 5 parameter PID type controller with delay
            mod.K = -mod.C*Gnm;
        case 5; % 4 parameter PID type controller with delay
            mod.K = -mod.C*Gnm*delay;
        case 6; % 4 parameter PID type controller with delay
            mod.K = -mod.C*Gnm*delay;
        case 7; % 4 parameter PID type controller with delay
            mod.K = -mod.C*Gnm*delay;
        otherwise
        disp('Non existing model number input');
    end

    % Calculate closed loop system responses
    me = [];
    mod.z = mod.G.zw + mod.G.zu*((eye(1)-mod.K*mod.G.yu)\mod.K*mod.G.yw);
    mod.y = mod.G.yw + mod.G.yu*((eye(1)-mod.K*mod.G.yu)\mod.K*mod.G.yw);
    try mod.y =  minreal(mod.y); catch me; end
    try mod.z =  minreal(mod.z); catch me; end
    disp(me);
    
     