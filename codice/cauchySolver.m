function b = cauchySolver(x,y,R,S,b)

% risolve sistema con matrice Cauchy-like A*v=b
% b può anche avere più di una colonna
% INPUT:
% R : generatore colonne A
% S : generatore righe A
% x : vettore colonna
% y : vettore riga
% b : vettore (o matrice) membro destro del sistema
% OUTPUT
% b : soluzione sistema

% problema: 
% la funzione sembra (nella pratica) più lenta di "\"
% nonostante l'andamento asintotico quadratico

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

n = length(x);
u = zeros(1,n);

for k = 1:n
    l = R * S(:,k) ./ ([y(1:k-1).';x(k:n)] - y(k));

    [valmax,imax] = max(abs(l(k:n)));
    if valmax == 0
        fprintf('attenzione: matrice numericamente singolare')
    end
    if imax > 1
        imax = imax + k - 1;
        l([k,imax]) = l([imax,k]);
        x([k,imax]) = x([imax,k]);
        b([k,imax],:) = b([imax,k],:);
        R([k,imax],:) = R([imax,k],:);
    end

    pivot = l(k);
    l(k) = -1;

    u(k+1:n) = R(k,:) * S(:,k+1:n) ./ (x(k) - y(k+1:n));

    R_pivot = R(k,:)/pivot;
    R(k,:) = 0;
    R = R - l * R_pivot;

    b_pivot = b(k,:)/pivot;
    b(k,:) = 0;
    b = b - l * b_pivot;

    S_pivot = S(:,k)/pivot;
    S(:,k+1:n) = S(:,k+1:n) - S_pivot * u(k+1:n);
end

end