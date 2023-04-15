function [h, gamma, beta] = rationalConverter(p,q)

% decompone funzione razionale p/q in fratti semplici
% INPUT:
% p, q : vettori polinomi
% OUTPUT: 
% h : vettore polinomio
% gamma : vettore coefficienti
% beta : vettore radici di q
% tali che :
% p(z)/q(z) = h(z) + gamma(1)/(z-beta(1)) + ... + gamma(n)/(z-beta(n))

beta = roots(q);

if length(p)<length(q)
    h = 0;
else
    [h,p] = deconv(p,q);
end

n = length(beta);
aux = zeros(n,1);

for i = 1:n
    vaux = beta(i) - beta( [1:(i-1), (i+1):n] );
    aux(i) = prod(vaux);
end
gamma = polyval(p,beta)./aux;

end