function y = tMatVec2(fT,x)

% prodotto matrice vettore per matrice di Toeplitz
% senza calcolare la fft associata alla matrice 
% T = toeplitz(a,b), che è passata come parametro fT

n = length(x);
y_aux = ifft(fT.*fft([x;zeros(n,1)]));
y = y_aux(1:n);

end
