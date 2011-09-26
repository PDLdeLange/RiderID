function [g,tau] = markovparm2(y,u,m1,m2)
% function [g,tau] = markovparm(y,u,m)
%
% Estimates a finite number of Markov parameters g(tau) based on input u
% and output y. the model is effectively an m*2+1 order FIR model.
% Since the process enhances a inverse operation, the estimation 
% may be slow for large datasets. 
%
n = length(y);
tau = (m1:m2)';
U = zeros(n,length(tau));
Y = zeros(n,1);

for i = 1:n
    for j = m1:m2
        if (i-j >= 1) && (i-j <= n)
            U(i,j-m1+1) = u(i-j);
        else 
            U(i,j-m1+1) = 0;
        end 
    end
    Y(i,:) = y(i);
end
g = (U'*U)\U'*Y;
