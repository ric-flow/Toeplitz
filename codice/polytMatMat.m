function Y =  polytMatMat(a,b,p,X)

% prodotto polinomio di matrice Toeplitz matrice
% INPUT:
% a : prima colonna matrice Toeplitz T
% b : prima riga matrice Toeplitz T
% p : polinomio
% X : matrice
% OUTPUT:
% Y : p(T)*X

a = a(:);
b = b(:);
d = length(p)-1; % grado del polinomio
[n,~] = size(X);

if d == 0
    Y = p*X;
end

if d == 1
    aux = fft([a;0;b(n:-1:2)]);
    Y = p(1) * tMatMat2(aux,X) + p(2)*X;
end

if d >= 2
    aux = fft([a;0;b(n:-1:2)]);
    Y = p(1) * tMatMat2(aux,X) + p(2)*X;
    for k = 3:d+1
        Y = tMatMat2(aux,Y) + p(k)*X;
    end
end

end
