function [G,B] = genCompress(G,B)

% dato un generatore (G,B) per A, lo comprime in modo che generi una 
% matrice vicina ad A, ma con displacement rank più basso
% INPUT:
% [G,B] : generatore originale
% OUTPUT:
% [G,B] : generatore compresso

% [n,m] = size(G);
[QG,RG] = qr(G,0);
[QB,RB] = qr(B,0);

aux = RG*RB';
[U,S,V] = svd(aux);

% d = diag(S);
% tol = d(1) * tol / n;
% pos = find(d<tol,1) - 1;
% if isempty(pos)
%     pos = m;
% end

pos = rank(S);

S1 = S(1:pos,1:pos);
U1 = U(:,1:pos);
V1 = V(:,1:pos);

G = QG * U1;
B = QB * V1 * S1;

end
