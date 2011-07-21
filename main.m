%% Davis data analysis and identification
clear all; close all; clc; addpath('functions','models');

% Loading data
try
    load('data/data.mat'); 
catch me
    disp(me);
    disp('No Data Detected -> dataprocessing...');
    raw = davisdat(); % Load davis Data
    bike = davisbike(raw.v); % Bicycle model from Davis
    fil = davisfilter(raw,0.1,[1 11600]);% Filtering davis data
    fir = firestimation(fil,2^9,0.2); % Finite impulse response
    sys = parametricmod(fir,bike,1); % Parametric model optimization
    res = results(fil,fir); % Results
    save('data/data.mat')
end

%% Plotting
close all;

fig(1) = firfig(fir,sys,1,1);
fig(2) = resfig(fil,res,2,1,[17 38]);