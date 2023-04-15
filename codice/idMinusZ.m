function G1 = idMinusZ(G)

% calcola (I-Z)*G

[n,m] = size(G);
G1 = G - [zeros(1,m); G(1:n-1,:)];

end