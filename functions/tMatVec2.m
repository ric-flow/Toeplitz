function y = tMatVec2(fT,x)

% prodotto matrice di Toeplitz vettore
% INPUT:
% fT : trasformata di Fourier associata alla matrice Toeplitz T
% x : vettore
% OUTPUT:
% y = T*x

n = length(x);
y_aux = ifft(fT.*fft([x;zeros(n,1)]));
y = y_aux(1:n);

end
