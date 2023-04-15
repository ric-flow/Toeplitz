function sol = T2CSolver(a,b,rhs)

% risolve sistema T*x=rhs con T Toeplitz convertendo in sistema A*y = w 
% con matrice A Cauchy-like
% INPUT:
% a : prima colonna di T
% b : prima riga di T
% rhs : vettore (o matrice) membro destro del sistema
% OUTPUT:
% sol : soluzione del sistema

a = a(:);
b = b(:);
if a(1)~=b(1)
    fprintf('attenzione, conflitto entrata (1,1) della matrice, la colonna vince il conflitto \n')
    b(1) = a(1);
end

n = length(a);

R = zeros(n,2);
S = zeros(n,2);

R(1,1) = 1;
S(n,2) = 1;

R(1,2) = 2 * a(1);
R(2:n,2) = a(2:n) + b(n:-1:2);
S(1:n-1,1) = a(n:-1:2) - b(2:n);

aux = unityRoots(2*n);
theta = aux(1:n);

R = ifft(R);
S = ifft(theta .* S);

x = aux(1:2:end);
y = aux(2:2:end);

rhs = ifft(rhs);
sol = cauchySolver(x,y,R,S,rhs);

sol = conj(theta) .* fft(sol) / n;
% sol = real(sol);

end