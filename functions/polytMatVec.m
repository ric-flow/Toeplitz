function y =  polytMatVec(a,b,p,x)

% prodotto polinomio di matrice Toeplitz vettore
% INPUT:
% a : prima colonna matrice Toeplitz T
% b : prima riga matrice Toeplitz T
% p : polinomio
% x : vettore
% OUTPUT:
% y : p(T)*x

a = a(:);
b = b(:);
x = x(:);
d = length(p)-1; % grado del polinomio
n = length(x);

if d == 0
    y = p*x;
end

if d == 1
    aux = fft([a;0;b(n:-1:2)]);
    y = p(1) * tMatVec2(aux,x) + p(2)*x;
end

if d >= 2
    aux = fft([a;0;b(n:-1:2)]);
    y = p(1) * tMatVec2(aux,x) + p(2)*x;
    for k = 3:d+1
        y = tMatVec2(aux,y) + p(k)*x;
    end
end

end
