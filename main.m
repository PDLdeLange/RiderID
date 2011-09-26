% Davis data analysis and identification
clear all; close all; clc; addpath('functions');

% Loading data
try % loading data file
    load('data/data.mat'); 
catch me
    %If data is not found then:
    disp(me);
    disp('No Data Detected -> Processing...');
    % Calculations:
    raw = davisdat('00105.mat'); % Load davis Data
    dat = prefilter(raw); % Filtering davis data
    set = settings(dat); % Measurement settings, etc.
    % Parametric and non parametric modelling
    npm = nonparametricmod(dat,set); % Nonparameteric modelling
    mod = parametricmod(npm,dat,set); % Parameteric modelling
    for i = 1:length(mod);
    disp(mod(i).vaf);
    disp(mod(i).X);
    fig(7) = covfig(mod(i),set,6);
    pause
    end

    % Perform Sensitivity Analysis for model i
    sym = symsens(mod(i));
    sen = sensitivity(sym,npm,mod(i));
    
    % Save data
    save('data/data.mat');
end


%% Dataplot

T = [18,43]; % Select nice datarange to compare
F = [0.1,10]; % Select nice datarange to compare

fig(1) = datfig(dat,2,T);
fig(2) = firfig(npm,mod(i),1,1);
fig(3) = firfig(npm,mod(i),2,2);
fig(4) = resfig(npm,dat,3,1,T);
fig(5) = resfig(npm,dat,4,2,T);
fig(6) = modfig(npm,dat,mod(i),6,1,T);
fig(7) = modfig(npm,dat,mod(i),7,2,T);

fig(7) = covfig(mod(i),set,6);
fig(8) = frffig(mod(i),npm,8,2,F);

%% Model comparison

bike = davisbike(set.v); % Bicycle model from Davis

G1 = mod(1);
G2 = mod(2);

G3 = optimal(bike);
G1_frf = squeeze(freqresp(G1.y,2*pi*npm.f)).';
G2_frf = squeeze(freqresp(G2.y,2*pi*npm.f)).';

figure(1);
bode(G1.y,G2.y,G3.y);

j = 2;
figure(8); clf;
leg = {'\phi','\delta'};
subplot(2,1,1); %limy = 0.02*[0 1];
    h = loglog(npm.f,abs([G1_frf(:,j) G2_frf(:,j)])); hold on;
    title(['Closed loop response: ' leg{j} '/w']); ylabel(['|' leg{j} '/w|']);
%     semilogx([sen.f_bw;sen.f_bw],repmat(limy',1,length(sen.f_bw)),':k'); hold on;
    legend(h,'PID','LQR');   xlim(F);% ylim(limy);
subplot(2,1,2);  %limy = 0.08*[0 1];
    semilogx(npm.f,unwrap(angle([G1_frf(:,j) G2_frf(:,j)]))); hold on;
%     title('Normalized parameter sensitivities');% ylabel('|d\phi/dXn|');
%     semilogx([sen.f_bw;sen.f_bw],repmat(limy',1,length(sen.f_bw)),':k'); hold on;
   ylabel(['angle(' leg{j} '/w)']);    xlim(F);% ylim(limy); 
xlabel('frequency (Hz)');
sdf('Latex');

%% Multisine results
Fs = 200;
Ts = 1/Fs;
P = 2^10;
D = 2^4;
N = P*D;
t = (1-1:N-1)'*Ts;
f = linspace(-1,1,P)'*Fs/2;
w_max = 40;
dF = 1/P/Ts;
f1 = 0.2;
f2 = 4.0;
[w,freqs] = idinput([P 1 D],'sine',...
[f1 f2]/(Fs/2),w_max*[-1 1],[2^3, 20, 1]);
disp(freqs*Fs/2/2/pi);

maxroll = 0.1745;
close all;
sim = modsim(npm,w,t);
y1max = max(abs(sim.y(:,1))); disp(y1max);
% plot(dat.y(:,1),'Color',0.4*ones(3,1)); hold on;
plot(sim.yhat(:,1),'Color',0.8*ones(3,1)); hold on;
plot(sim.y(:,1),'Color',0.4*ones(3,1)); hold on;


%% Schommeldesign

parms.g =    9.81;
parms.c =   10.00;
parms.m =    5.00;
parms.e =    0.30;
parms.h =    0.75;
parms.d =    0.65;
sys  = schommel(parms);
sys  = slider(parms);
disp(eig(sys.theta));

theta = lsim(sys.theta,w,t);
plot(t,theta);
disp(max(abs(theta)));



%% Sensitivity plot
figure(8); clf;
subplot(2,1,1); limy = 0.02*[0 1];
    h = semilogx(sen.f,abs(squeeze(sen.dydXn(1,:,:)))); hold on;
    title('Normalized parameter sensitivities'); ylabel('|d\phi/dXn|');
    semilogx([sen.f_bw;sen.f_bw],repmat(limy',1,length(sen.f_bw)),':k'); hold on;
    legend(h,'X1','X2','X3','X4','X5','X6','X7'); ylim(limy);
subplot(2,1,2);  limy = 0.08*[0 1];
    semilogx(sen.f,abs(squeeze(sen.dydXn(2,:,:)))); hold on;
    title('Normalized parameter sensitivities'); ylabel('|d\phi/dXn|');
    semilogx([sen.f_bw;sen.f_bw],repmat(limy',1,length(sen.f_bw)),':k'); hold on;
    ylabel('angle(d\delta/dXn)'); ylim(limy);
xlabel('frequency (Hz)');
sdf('Latex');
for k = 1:2
    disp((sen.cov_bw(:,:,k)))
end
    
%% Disturbance Spectrum Svv
figure(9); clf;
h = loglog(npm.f,npm.Svv); hold on;
legend(h,'\phi','\delta');
title('Disturbance Spectrum (Svv)'); hold on; box on;
xlabel('frequency (Hz)');
ylabel('|Svv|');
sdf('Latex');

%% Latex output
X0table = cell2mat({mod.X0}'); 
Xtable = cell2mat({mod.X}'); disp(Xtable);
vaf = cell2mat({mod.vaf}); disp(vaf);
s = mlatex([Xtable vaf'], '%.2f');

