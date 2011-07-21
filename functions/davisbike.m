function sys = davisbike(v)
% sys = davisbike(v)
%
% input v: forward velocity of the bicycle.
% output sys: state space model of the bicycle.
%

    % Load parameters into a structure
    addpath('functions','models');
    par = par_text_to_struct('data/bicycleparms.m');
    % Generates state space matrices
    [aMat, bMat, ~, ~] = whipple_pull_force(par, v);
    
    i = [9;11;4;7]; % Select state: phid, deltad, phi and delta.
    
    % Set up State Space matrices
    A = aMat(i,i);
    B = bMat(i,:);
    C = eye(4);
    D = 0;
    
    % Combine A,B,C and D matrices into a state space object.
    sys = ss(A,B,C,D);
    
end
