function Y = rationalEval(X,h,gamma,beta)

% h, gamma, beta descrivono la funzione razionale decomposta
% r(z) = h(z) + gamma(1)/(z-beta(1)) + ... + gamma(s)/(z-beta(s))
% vengono passati da rationalconverter (vedere funzione)
% la funzione restituisce Y = r(X)

Y = polyvalm(h,X);
n = size(X,1);
s = length(beta);

for i = 1:s
    Y = Y + gamma(i) * inv( X - beta(i) * eye(n) );
end

end
