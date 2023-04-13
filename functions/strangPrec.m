function c = StrangPrec2(a,b)

% data la prima colonna e la prima riga di una matrice di Toeplitz
% genera la prima colonna del precondizionatore circolante di Strang

a = a(:);
b = b(:);
n = length(a);

m = floor(n/2);
if n == 2*m
    c = [a(1:(m)); (a(m+1)+b(m+1))/2; b(m:-1:2)];
else
    c = [a(1:(m+1)); b((m+1):-1:2)];
end
end