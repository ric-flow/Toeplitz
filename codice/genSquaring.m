function [G2,B2] = genSquaring(G,B)

% INPUT:
% G,B generatore della matrice A
% OUTPUT:
% G2,B2 generatore della matrice A^2

[m,n]=size(G);
if (n>m)
  G=G';
end
[m,n]=size(B);
if (n>m)
  B=B';
end

[n,m] = size(G);
G_sum = cumsum(G);
B_sum = cumsum(B);

G2 = zeros(n, 2*m+1);
B2 = zeros(n, 2*m+1);
G2_sum = zeros(n,m+1);
B2_sum = zeros(n,m+1);

for k = 1:m
    G2_sum(:,k) = genMatVec(G,B,G_sum(:,k));
end
G2(:,m+1:2*m) = G;
G2_sum(:,m+1) = -genMatVec(G,B,ones(n,1));

B2(:,1:m) = B;
for k = 1:m
    B2_sum(:,k) = genMatVec(B,G,B_sum(:,k));
end
B2_sum(:,m+1) = genMatVec(B,G,ones(n,1));

G2(:,[1:m, end]) = idMinusZ(G2_sum);
B2(:,m+1:end) = idMinusZ(B2_sum);

end