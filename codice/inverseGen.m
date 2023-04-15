function [G,B] = inverseGen(a,b,toepSolver)

% calcola generatore (Stein) dell'inversa di T Toeplitz
% INPUT:
% a : prima colonna di T
% b : prima riga di T
% toepSolver : handle function tale che
% x = toepSolver(a,b,y) risolve il sistema T*x = y
% OUTPUT:
% [G,B] generatore di inv(T)

a = a(:);
b = b(:);

n = length(a);
e = ones(n,1);

a_aux = cumsum(a);
b_aux = cumsum([0; b(2:n)]);

G(:,1) = e - toepSolver(a,b,a_aux) ;
G(:,2) = toepSolver(a,b,e) ;
B(:,1) = conj(G(n:-1:1,2)) ;
B(:,2) = e - toepSolver(b',a',b_aux) ;

G = idMinusZ(G);
B = idMinusZ(B);

end
