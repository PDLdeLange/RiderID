function Gs = lti2sym(G)
% function Gsym = lti2sym(G)
%
% Converts MIMO LTI models to symbolic transfer function expresions, where
% 's' denotes the Laplace argument. This function is not tested with 
% time delays. 
%
% Peter de Lange 2011 (pdldelange@gmail.com)

syms s;

[ni,nj] = size(G);
Gs = sym(zeros(ni,nj));
for i = 1:ni
    for j = 1:nj
    [num,den] = tfdata(G(i,j));
    Gs(i,j) = poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s);
    end
end
