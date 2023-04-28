function [G,B] = simpleRatGen(a,b,h,gamma,beta,toepSolver)

% calcola generatore (Stein) di funzione razionale 
% (decomposta in fratti semplici) valutata in T Toeplitz
% INPUT:
% a : prima colonna di T
% b : prima riga di T
% h, gamma, beta : vettori descrivono r(z) decomposta
% r(z) = h(z) + gamma(1)/(z-beta(1)) + ... + gamma(s)/(z-beta(s))
% toepSolver : handle function tale che
% x = toepSolver(a,b,y) risolve il sistema T*x = y
% OUTPUT:
% [G,B] generatore di r(T)

a = a(:);
b = b(:);
n = length(a);

if h ~= 0
    [Gh,Bh] = polynomialGen(a,b,h,0);
end

% generatore fratti semplici ottimizzato sfruttando che p,q sono reali
s = length(beta);
G_rat = zeros(n,2*s);
B_rat = zeros(n,2*s);
ibeta = imag(beta);
index_real = find(ibeta==0)';
for i = index_real
    ai = a;
    ai(1) = a(1) - beta(i);
    b(1) = ai(1);
    [G_rat(:,(2*i-1):(2*i)),B_rat(:,(2*i-1):(2*i))] = inverseGen(ai,b,toepSolver);
    G_rat(:,(2*i-1):(2*i)) = gamma(i)*G_rat(:,(2*i-1):(2*i));
end
index_notreal = find(ibeta ~= 0)';
for i = index_notreal(1:2:end)
    ai = a;
    ai(1) = a(1) - beta(i);
    b(1) = ai(1);
    [G_rat(:,(2*i-1):(2*i)),B_rat(:,(2*i-1):(2*i))] = inverseGen(ai,b,toepSolver);
    G_rat(:,(2*i-1):(2*i)) = gamma(i)*G_rat(:,(2*i-1):(2*i));
end
for i = index_notreal(2:2:end)
    G_rat(:,(2*i-1):(2*i)) =  conj(G_rat(:,(2*i-3):(2*i-2)));
    B_rat(:,(2*i-1):(2*i)) =  conj(B_rat(:,(2*i-3):(2*i-2)));
end

if h ~= 0
    G = [Gh,G_rat];
    B = [Bh,B_rat];
else
    G = G_rat;
    B = B_rat;
end

end