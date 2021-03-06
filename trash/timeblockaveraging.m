clear all; close all; clc; 

dav = davisdat();

% Davis results check
figure(1);
subplot(3,1,1);
    plot(dav.t,dav.y(:,1),'k'); hold on;
    ylabel('Steering Angle (rad)');
    title('Raw measurements');
subplot(3,1,2);
    plot(dav.t,dav.y(:,2),'k'); hold on;
    ylabel('Roll Angle (rad)');
subplot(3,1,3);
    plot(dav.t,dav.w,'k');
    ylabel('Force input (N)');
    xlabel('Time (s)');
sdf('Latex');


%% Filtering

% Divide the measurements into D blocks  of P samples.
P = 2^10; shift = 10; D = 7;

% Aprropriate time and frequency vector
t = linspace(0,(P-1)/dav.Fs,P)';
f = (0:P/2-1)'/t(end);

% Simple but ugly algorithm for selecting blocks based on treshold on
% force input. 
u = zeros(P,D); y1 = zeros(P,D); y2 = zeros(P,D);
j = 1; i = 1;
while i < dav.N; % For every sample
    if abs(dav.w(i)) > 20; % Check wether input treshhold is exceeded
        sel = (1:P)+i-shift;
        w(:,j) = dav.w(sel,1); % Force selection
        y1(:,j) = dav.y(sel,1); % Roll angle selection
        y2(:,j) = dav.y(sel,2); % Steering angle selection
        disp([i,j])
        i = i + P; % Skip first P samples so that the blocks don't overlap.
        j = j + 1; % Next block
    end
    i = i+1; % Next sample
    if j>D; break; end; % Enough is enough!!!!!
end

%% Identification toolbox





% 
% data = iddata(y,w,dav.Fs^-1);
% 
% m = armax(data(:,1),[4 4 4 4]);
% m = pem(data)
% ymod = sim(m,w);


%% Plotting

% Datacheck
figure(2);
subplot(3,1,1);
title('Filtered Measurments');
    plot(t,y1,':'); hold on;
    plot(t,mean(y1,2),'k-');
%     plot(t,mean(ymod,2),'k-');
    ylabel('Roll angle (rad)');
subplot(3,1,2);
    plot(t,y2,':'); hold on;
    plot(t,mean(y2,2),'k-');
    ylabel('Steering Angle (rad)');
subplot(3,1,3);
    plot(t,u,':'); hold on;
    plot(t,mean(u,2),'k-');
    ylabel('Force (N)');
xlabel('Time (s)');
sdf('Latex');



%% SYSID



par = par_text_to_struct('bicycleparms.m');

% State space equations of the bicycle.
[Atmp, Btmp, ~, ~] = whipple_pull_force(par,dav.v);
T = zeros(4,11); T(1,9) = 1; T(2,11) = 1; T(3,4) = 1; T(4,7) = 1;
A = T*Atmp*T'; B = T*Btmp; C = eye(4); clear Atmp Btmp T;
bike = ss(A,B,C,0); clear A B C D;

%% Fitting

X0 = ones(1,8);
% Xn =[  -11.9903   1.8494  -39.3832   23.2078];
% Xn =[   -0.0026   12.5181  -12.3760];
% Xn =[ -10 -11.9903  -39.3832 ];
Xn =[0 0  -11.9903   1.8494  -39.3832   23.2078 0 0];
Xn =[-0.3950   -0.1893    0.0000   -2.5450    4.7676    1.9137   -1.2782    0.6145];
efunc = @(Xn)errorfunc(Xn,X0,w,y,t,bike);

% Ky      =   .1021;
% Kpsi    =   .1651; 
% Kphi    =  8.5460; 
% Kphid   = -0.0750;
% Kdelta  = 60.0000;
% Xn = [Kphid 0 Kphi Kdelta];

Xn = lsqnonlin(efunc,Xn);

sys = riderfunc(Xn,X0,bike);
ymod = lsim(sys.y,w,t);
figure(5); cla;
plot(t,y,'-'); hold on;
plot(t,ymod,':');
disp(eig(sys.y));

%% TRASH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% par = par_text_to_struct('bicycleparms.m');
% 
% % State space equations of the bicycle.
% [Atmp, Btmp, ~, ~] = whipple_pull_force(par,dav.v);
% T = zeros(4,11); T(1,9) = 1; T(2,11) = 1; T(3,4) = 1; T(4,7) = 1;
% A = T*Atmp*T'; B = T*Btmp; C = eye(4); clear Atmp Btmp T;
% bike = ss(A,B,C,0); clear A B C D;
% % 
% %% Bicycle/rider system
% 
% [sys.G] = bike(:,1:2);
% % [sys] = optimal(sys);
% 
% %% Davis Datasim
% 
% % M = inv(sys.G.b(1:2,1:2)); % Recover Mass Matrix
% % f = bike.b(1:2,3);
% % u(:,1) = M (1,:)*f(:)*dav.force; u(:,2) =M (2,:)*f(:)*dav.force;

% %% Frequency domain methods
% 
% wu = zeros(P,1); wu(1:(P/2^3),:) = 1;
% wy = exp(log(0.1)*((1:P)-1)/((P-1)))';% decaying exponential window
% 
% U  = (fft(repmat(wu,1,D).*u));  U  = U(1:P/2,:);
% Y3 = (fft(repmat(wy,1,D).*y1)); Y3 = Y3(1:P/2,:);
% Y4 = (fft(repmat(wy,1,D).*y2)); Y4 = Y4(1:P/2,:);
% 
% Sy1y1 = 1/D*(conj(Y3).*Y3);
% Sy2y2 = 1/D*(conj(Y4).*Y4);
% Suy1 = 1/D*(conj(U).*Y3);
% Suy2 = 1/D*(conj(U).*Y4);
% Suu = 1/D*(conj(U).*U);
% 
% G3 = mean(Suy1./Suu,2);
% G4 = mean(Suy1./Suu,2);
% 
% %% FOUT:  Gc ipv. bike:  Gmod = squeeze(freqresp(ss(bike.a,bike.b(:,3),bike.c(3,:),0),2*pi*f));
% 
% coh3 = mean(abs(Suy1),2).^2./(mean(Suu,2).*mean(Sy1y1,2));
% coh4 = mean(abs(Suy2),2).^2./(mean(Suu,2).*mean(Sy2y2,2));
% 
% % Plotting
% figure(3);
% subplot(3,1,1);
% loglog(f,abs(G3),'k'); hold on;
% loglog(f,abs(Gmod),'r:');
% subplot(3,1,2);
% semilogx(f,phase(G3),'k'); hold on;
% loglog(f,phase(Gmod),'r:');
% ylim(2*pi*[-1 1]);
% subplot(3,1,3);
% semilogx(f,coh3); ylim([0 1]);
% sdf('Latex');
% 
% 
% % e=exp(-time*a);		% decaying exponential window
% % e1=1./e;		% rising exponential window
% 
% 
