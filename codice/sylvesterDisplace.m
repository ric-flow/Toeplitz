function DeltaX = sylvesterDisplace(X)

% displacement di tipo Sylvester applicato alla matrice X
% INPUT:
% X : matrice
% OUTPUT:
% DeltaX = Z1*X - X*Zminus1

DeltaX = [X(end,:); X(1:end-1,:)] - [X(:,2:end), -X(:,1)];

end
