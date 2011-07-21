function dav = davisdat()
% Creates structure with some of the Davis measurement Data. 

dat = load('data/00105.mat');

Fs = dat.NISampleRate; Fs = double(Fs);

v = dat.ForwardSpeed; 

y(:,1) = dat.RollAngle;
y(:,2) = dat.SteerRate;
y(:,3) = dat.RollAngle;
y(:,4) = dat.SteerAngle;
f      = dat.PullForce;

N = size(y,1);

t = (0:N-1)'./Fs;
% Output

dav.Fs = Fs;
dav.v = mean(v);
dav.t = t;
dav.y = y;
dav.N = N;
dav.force = f;

