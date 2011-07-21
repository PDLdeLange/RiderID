function [H] = bikesys(v,datapath)
% function sys = bikesys{v,datapath)
%
% Outputs a state space model of the bicycle equations
%
% Inputs:
% - v; forward speed [double]
% - g; gravitational acceleration [double]
% - datapath; path to bicycle .mat file (string)
%
% Outputs:
% - sys; state space respresentation of the bicycle equations (ss-object)
    
    load(datapath);
    M = roundn(M0,-12);
    H = struct([]);
    g = 9.81;
    for i = 1:length(v);
    
        
        C = roundn((C1*v(i)),-12);
        K = roundn((K0*g + K2*v(i)^2),-12);

        sys.A = [-M\C -M\K; eye(3) zeros(3);];
        sys.B = [ M\eye(3); zeros(3)];
        sys.C = eye(6);
        sys.D = zeros(6,3);

        H(i).v = v(i);
        H(i).G = minreal(ss(sys.A,sys.B,sys.C,sys.D));
        
    end

  