function [G,B] = rationalGen(a,b,p,q,toepSolver)

% va sistemato, deve esserci un errore da qualche parte
% per ora utilizzare simpleRatGen

a = a(:);
b = b(:);
n = length(a);

dp = length(p) - 1;
dq = length(q) - 1;

G_aux = zeros(2*n,2*(dp+dq)+1);
B_aux = zeros(2*n,2*(dp+dq)+1);

[G_aux(1:n,1:2*dq),B_aux(1:n,1:2*dq)] = polynomialGen(a,b,q,1);
[G_aux(1:n,2*dq+1:2*(dq+dp)),B_aux(n+1:2*n,2*dq+1:2*(dq+dp))] = polynomialGen(a,b,p,1);
G_aux(n+1:2*n,end) = - ones(n,1);
B_aux(1:n,end) = ones(n,1);

[G_aux,B_aux] = genCompress(G_aux,B_aux);

G = G_aux(n+1:2*n,:) + toepSolver(a,b,q,G_aux(1:n,:));
B = B_aux(n+1:2*n,:) - polytMatMat(b,a,p, toepSolver(b,a,q,B_aux(1:n,:)) );

G = idMinusZ(G);
B = idMinusZ(B);

end
