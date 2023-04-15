function zeta = unityRoots(n)

% Restituisce il vettore (colonna) delle radici complesse n-esime di 1
% INPUT:
% n
% OUTPUT:
% zeta : vettore radici n-esime

zeta = exp( complex( 0, 2 * pi * (0:n-1).' / n ) );

end