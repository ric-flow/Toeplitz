function Y =  polytMatMat2(fT,p,X)

% prodotto polinomio di matrice Toeplitz matrice
% INPUT:
% fT : trasformata di Fourier associata alla matrice Toeplitz T
% p : polinomio
% X : matrice
% OUTPUT:
% Y : p(T)*X

d = length(p)-1; % grado del polinomio

if d == 0
    Y = p*X;
end

if d == 1
    Y = p(1) * tMatMat2(fT,X) + p(2)*X;
end

if d >= 2
    Y = p(1) * tMatMat2(fT,X) + p(2)*X;
    for k = 3:d+1
        Y = tMatMat2(fT,Y) + p(k)*X;
    end
end

end
