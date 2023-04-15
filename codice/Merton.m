function [a, b, T] = Merton(n)

% crea matrice di Toeplitz associata alla discretizzazione della 
% PIDE (Partial Integro-Differential Equation) del Merton Model
% INPUT:
% n : taglia desiderata per la matrice
% OUTPUT:
% a : prima colonna di T
% b : prima riga di T
% T : matrice del Merton Model

% inizializzazione parametri
xi_min = -2;
xi_max = 2;
Delta_xi = (xi_max - xi_min) / (n+1);
nu = 0.25;
r = 0.05;
lambda = 0.1;
mu = -0.9;
sigma = 0.45;

kappa = exp( mu + sigma^2/2) - 1;

theta = @(eta) exp( - (eta-mu).^2 / (2*sigma^2) ) / (sqrt(2 * pi) * sigma);

% Matrice tridiagonale corrispondente alla parte differenziale della PIDE
mid = nu^2/Delta_xi^2;
off = (2*r - 2*lambda*kappa - nu^2)/(4*Delta_xi);
D = diag((mid/2 - off) * ones(n-1,1), -1 ) + ...
    diag( (-mid - r - lambda) * ones(n,1) ) + ...
    diag((mid/2 + off) * ones(n-1,1), 1);

% Matrice di Toeplitz corrispondente alla parte integrale della PIDE
Ic = theta( ( 0:-1:(1-n) )' * Delta_xi );
Ir = theta( ( 0:(n-1) )  * Delta_xi );
I = Delta_xi * toeplitz(Ic, Ir);

% Assembla la matrice
T = D + lambda * I;
a = T(:,1);
b = T(1,:).';

end