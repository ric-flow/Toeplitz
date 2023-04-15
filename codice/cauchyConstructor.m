function A = cauchyConstructor(R,S,x,y)

% costruisce matrice Cauchy-like
% INPUT:
% R : generatore colonne
% S : generatore righe
% x : vettore colonna
% y : vettore riga
% OUTPUT
% A : matrice Cauchy-like che rispetta il displacement
% D_x * A - A * D_y = R * S'

[m,n]=size(R);
if (n>m)
  R=R.';
end

[m,n]=size(S);
if (m>n)
  S=S';
end

[m,n]=size(x);
if (n>m)
  x=x.';
end

[m,n]=size(y);
if (m>n)
  y=y.';
end

A=R*S./(x-y);

end