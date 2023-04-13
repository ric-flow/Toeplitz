function Y = tMatMat2(fT,X)

% prodotto matrice di Toeplitz matrice
% INPUT:
% fT : trasformata di Fourier associata alla matrice Toeplitz
% X : matrice
% OUTPUT:
% Y = T*X prodotto matrice matrice

[n,m] = size(X);
Y_aux = ifft(fT.*fft([X;zeros(n,m)]));
Y = Y_aux(1:n,:);

end
