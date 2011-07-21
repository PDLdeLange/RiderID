function [sys] = riderfunc(X,bike)

% Bike
sys.Gyu =  bike(3:4,2);
sys.Gyw =  bike(3:4,3);
sys.Gzu = -bike(3,2);
sys.Gzw = -bike(3,3);

% State gains
sys.gains = -X;
Gnm = tf(900,[1 2*.707*30  900]);
s = tf('s');
sys.K = sys.gains*[eye(2)/s; eye(2)*s; eye(2)]*Gnm*exp(-0.00*s);

sys.z = minreal(sys.Gzw + sys.Gzu*((eye(1)-sys.K*sys.Gyu)\sys.K*sys.Gyw));
sys.y = minreal(sys.Gyw + sys.Gyu*((eye(1)-sys.K*sys.Gyu)\sys.K*sys.Gyw));