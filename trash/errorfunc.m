function e = errorfunc(Xn,X0,u,y,t,bike)

[sys] = riderfunc(Xn,X0,bike);

ymod = lsim(sys.y,u,t);

e = y(:,2)-ymod(:,2);