function [G,B] = monomialGen(a,b,k)

% calcola un generatore Stein per T^k
% INPUT:
% a : prima colonna di T
% b : prima riga di T
% k : grado del monomio T^k
% OUTPUT:
% [G,B] : generatore per T^k

a = a(:);
b = b(:);
n = length(a);

if k == 0
    G = eye(n,1);
    B = eye(n,1);
    return
end

if k == 1
    G = zeros(n,2);
    G(:,1) = a;
    G(:,2) = eye(n,1);
    B(:,1) = eye(n,1);
    B(:,2) = [0; b(2:n)];
    return
end

c = [a; 0; b(n:-1:2)];
c = fft(c);

G = zeros(n,2*k);
B = zeros(n,2*k);

% inizializza variabili ausiliarie
a_aux = cumsum(a);
b_aux = cumsum([0; b(2:n)]);
b_aux = b_aux(n:-1:1);

% riempie G
G(:,1) = -b_aux;
G(:,2) = ones(n,1);
G(:,3) = tMatVec2(c,-b_aux);
G(:,4) = a_aux + b_aux;

for j = [5:(2*k-2), 2*k]
    G(:,j) = tMatVec2(c,G(:,j-2));
end

% la penultima colonna è diversa dalle altre, vedere teoria
G(:,2*k-1) = G(:,2*k-3) + G(:,2*k);
G(:,2*k-1) = tMatVec2(c,G(:,2*k-1));

% ricalcola variabili ausiliarie
b_aux = b_aux(n:-1:1);
a_aux = a_aux(n:-1:1);

% calcola la fft del polinomio associato a T^*
c = [b; 0; a(n:-1:2)];
c = fft(c);

% riempie B
B(:,2*k) = b_aux;
B(:,2*k-1) = ones(n,1);
B(:,2*k-2) = tMatVec2(c,b_aux);
B(:,2*k-3) = b_aux + a_aux;

for j = (2*k-4):(-1):(1)
    B(:,j) = tMatVec2(c,B(:,j+2));
end

G = idMinusZ(G);
B = idMinusZ(B);

end
