function c = StrangPrec2(a,b)

% genera precondizionatore circolante di Strang
% INPUT:
% a : prima colonna matrice Toeplitz T
% b : prima riga matrice Toeplitz T
% OUTPUT:
% precondizionatore circolante di Strang C(T)

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