function y = norm1T(a,b)

% INPUT:
% a : prima colonna di T
% b : prima riga di T
% OUTPUT: 
% norma 1 della matrice T in tempo O(n)

n= length(a);
aux = sum(abs(a));
y = aux;

for i = 2:n
    aux = aux - abs(a(n-i+2)) + abs(b(i));
    if aux > y
        y = aux;
    end
end

end
