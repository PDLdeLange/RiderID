%% Davis data analysis and identification
clear all; close all; clc; addpath('functions');

% Loading data
try % loading data file
    load('data/data.mat'); 
catch me % If data is not found then do ...
    disp(me);
    disp('No Data Detected -> Processing...');
    raw = davisdat(); % Load davis Data
    fil = davisfilter(raw,0.1);% Filtering davis data
    fir = firestimation(fil,2^9,0.1); % Finite impulse response
    etf = frfestimation(fil);
    arx = arxestimation(fil);
    res = results(fil,fir); % Results
    bike = davisbike(raw.v); % Bicycle model from Davis
    mod = parametricmod(fir,bike,1); % Parametric model optimization
    for i = 1:length(mod); % for each model
        disp([num2str(i) '/' num2str(length(mod))])
        [frf(i),sym(i)] = sensitivity(mod(i)); %#ok<*SAGROW>
    end
    save('data/data.mat');
end

%% Plotting
close all;

i = 6;
% System i
fig(1) = firfig(fir,mod(i),1,1);
fig(2) = firfig(fir,mod(i),2,2);
fig(3) = resfig(fil,res,3,1,[res.t(1), res.t(end)]);
% fig(4) = resfig(fil,res,4,2,[res.t(1), res.t(end)]);
fig(5) = frffig(etf,mod(i),5);
fig(6) = covfig(mod(i),6);
fig(7) = sensfig(frf(i),7);

% Results
Xtable = cell2mat({mod.X}'); disp(Xtable);
vaf = cell2mat({mod.vaf}); disp(vaf);



