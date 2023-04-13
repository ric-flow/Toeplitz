function Y = tMatMat(a,b,X)

% prodotto matrice di Toeplitz matrice
% INPUT: 
% a : prima colonna matrice Toeplitz T
% b : prima riga matrice Toeplitz T
% X : matrice
% OUTPUT:
% Y : T*X

a = a(:);
b = b(:);
[n,m] = size(X);

aux = [a;0;b(n:-1:2)];
Y_aux = ifft(fft(aux).*fft([X;zeros(n,m)]));
Y = Y_aux(1:n,:);

end
