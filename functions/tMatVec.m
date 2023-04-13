function y = tMatVec(a,b,x)

% prodotto matrice di Toeplitz vettore
% INPUT: 
% a : prima colonna matrice Toeplitz T
% b : prima riga matrice Toeplitz T
% x : vettore
% OUTPUT:
% y : T*x

a = a(:);
b = b(:);
x = x(:);
n = length(x);

aux = [a;0;b(n:-1:2)];
y_aux = ifft(fft(aux).*fft([x;zeros(n,1)]));
y = y_aux(1:n);

end