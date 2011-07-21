function [g,tau] = markovparm(y,u,m)
% function [g,tau] = markovparm(y,u,m)
%
% Estimates a finite number of Markov parameters g(tau) based on input u
% and output y. the model is effectively an m*2+1 order FIR model.
% Since the process enhances a inverse operation, the estimation 
% may be slow for large datasets. 
%
n = length(y);
tau = (-m:m)';
U = zeros(n-2*m,2*m+1);
Y = zeros(n-2*m,1);
for i = m+1:n-m
    for j = -m:m
        U(i-m,j+m+1) = u(i-j);
        Y(i-m,:) = y(i);
    end
end
g = U\Y;