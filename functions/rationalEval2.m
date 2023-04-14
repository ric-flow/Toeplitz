function w = rationalEval2(X,h,gamma,beta,v)

% h, gamma, beta descrivono la funzione razionale decomposta
% r(z) = h(z) + gamma(1)/(z-beta(1)) + ... + gamma(s)/(z-beta(s))
% vengono passati da rationalconverter (vedere funzione)
% la funzione restituisce w = r(X)*v

Y = polyvalm(h,X);
s = length(beta);
n = size(X,1);
w = Y*v;

for i = 1:s
    w = w + gamma(i) * (( X - beta(i) * eye(n) ) \ v);
end


end