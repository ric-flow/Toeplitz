function NablaX = steinDisplace(X)

% displacement di tipo Stein applicato alla matrice X
% INPUT:
% X : matrice
% OUTPUT:
% NablaX = X - Z*X*Z

X(2:end,2:end) = X(2:end,2:end) - X(1:end-1,1:end-1);
NablaX = X;

end
