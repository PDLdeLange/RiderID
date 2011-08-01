%% Davis data analysis and identification
clear all; close all; clc; addpath('functions','models');

% Loading data
try % loading data file
    load('data/data.mat'); 
catch me % If data is not found then do ...
    disp(me);
    disp('No Data Detected -> Processing...');
    raw = davisdat(); % Load davis Data
    fil = davisfilter(raw,0.1);% Filtering davis data
    fir = firestimation(fil,2^9,0.2); % Finite impulse response
     bike = davisbike(raw.v); % Bicycle model from Davis
    [sys,opt] = parametricmod(fir,bike,1); % Parametric model optimization
    res = results(fil,fir); % Results
    save('data/data.mat')
end

%% Plotting
close all;

fig(1) = firfig(fir,sys,1,1);
fig(2) = firfig(fir,sys,2,2);
fig(3) = resfig(fil,res,3,1,[res.t(1), res.t(end)]);[17 38];
fig(4) = resfig(fil,res,4,2,[res.t(1), res.t(end)]);[17 38];

