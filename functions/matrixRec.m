function X = matrixRec(NablaX)

% ricostruisce una matrice dato il suo displacement di tipo Stein
% INPUT:
% NablaX : displacement della matrice X
% OUTPUT:
% X

[n,m] = size(NablaX);
X = NablaX;

for k = 2:m
    X(:,k) = X(:,k) + [0; X(1:n-1,k-1)];
end

end