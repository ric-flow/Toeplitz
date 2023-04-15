function [A,dr] = expmT(a,b,p,q,theta,toepSolver)

% calcola l'esponenziale della matrice di Toeplitz T = toeplitz(a,b)
% tramite una tecnica di scaling e squaring strutturato, con approssimante
% razionale p/q
% INPUT: 
% a : prima colonna matrice Toeplitz T
% b : prima riga matrice Toeplitz T
% p,q : polinomi, in modo che p/q è approssimante scelto per l'esponenziale
% (di solito si scelgono approssimanti di Padé diagonali)
% theta : parametro di scaling (scala finché la norma è minore di theta)
% toepSolver : handle function tale che
% x = toepSolver(a,b,y) risolve il sistema T*x = y
% OUTPUT:
% A : matrice A = exp(T)
% dr : vettore dei displacement rank ai vari passi di squaring

n1 = norm1T(a,b);
s = 0;
while n1 > theta
    n1 = n1/2;
    s = s+1;
end
m = 2^s;
a = a/m;
b = b/m;

[h, gamma, beta] = rationalConverter(p,q);

[G,B] = simpleRatGen(a,b,h,gamma,beta,toepSolver);
[G,B] = genCompress(G,B);
G = real(G);
B = real(B);
dr = zeros(s+1,1);
dr(1) = size(G,2);

for k = 1:s
    [G,B] = genSquaring(G,B);
    [G,B] = genCompress(G,B);
    dr(k+1) = size(G,2);
end

A = matrixRec(G*B');

end