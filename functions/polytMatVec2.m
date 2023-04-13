function y =  polytMatVec2(fT,p,x)

% prodotto polinomio di matrice Toeplitz vettore
% INPUT:
% fT : trasformata di Fourier associata alla matrice Toeplitz T
% p : polinomio
% x : vettore
% OUTPUT:
% y : p(T)*x

x = x(:);
d = length(p)-1; % grado del polinomio

if d == 0
    y = p*x;
end

if d == 1
    y = p(1) * tMatVec2(fT,x) + p(2)*x;
end

if d >= 2
    y = p(1) * tMatVec2(fT,x) + p(2)*x;
    for k = 3:d+1
        y = tMatVec2(fT,y) + p(k)*x;
    end
end

end