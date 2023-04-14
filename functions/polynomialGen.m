function [G,B] = polynomialGen(a,b,p)

% calcola un generatore Stein per p(T)
% INPUT:
% a : prima colonna di T
% b : prima riga di T
% p : polinomio
% OUTPUT:
% [G,B] : generatore per p(T)

a = a(:);
b = b(:);
n = length(a);

k = length(p) - 1; % grado del polinomio

if k == 0
    G = zeros(n,1);
    G(1) = p;
    B = eye(n,1);
    return
end

[G_aux,B] = monomialGen(a,b,k);

G = p(1) * G_aux;

for j = 1 : k - 1
    G_aux(:,2*k-3) = G_aux(:,2*k-3) + G_aux(:,2*k);
    G_aux(:,2*j+1:2*k) = G_aux(:,2*j-1:2*k-2);
    G_aux(:,1:2*j) = zeros(n,2*j);
    G = G + p(j+1) * G_aux;
end

G(1,2*k-1) = G(1,2*k-1) + p(k+1);

end
