function y = genMatVec(G,B,x)

% calcola il prodotto y = A*x in tempo subquadratico,
% sfruttando un generatore
% INPUT:
% [G,B] : generatore di A
% x : vettore
% OUTPUT:
% y = A*x

[m,n]=size(G);
if (n>m)
  G=G';
end
[m,n]=size(B);
if (n>m)
  B=B';
end
[m,n]=size(x);
if (n>m)
  x=x.';
end

B = conj(B);

[n,m] = size(G);
B_aux = [B(1,:); zeros(n,m); B(n:-1:2,:)];
x_aux = [x; zeros(n,1)];
v_aux = ifft(fft(B_aux).*fft(x_aux));
v = v_aux(1:n,:);
G_aux = [G; zeros(n,m)];
v_aux = [v; zeros(n,m)];
w_aux = ifft(fft(G_aux).*fft(v_aux));
w = w_aux(1:n,:);
y = sum(w,2);

end
