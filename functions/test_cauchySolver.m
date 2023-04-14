% test per il solver di sistemi cauchy
% la nostra implementazione non è molto veloce
% c'è comunque un guadagno nei tempi di cpu per matrici grandi (4000 o più)

n = 4000;
dr = 4;
b_size = 7;
b = rand(n,b_size);
R = rand(n,dr);
S = rand(n,dr);
x = 1:2:2*n;
y = 2:2:2*n;
A = cauchyConstructor(R,S,x,y);

sol1 = cauchySolver(x,y,R,S,b);
sol2 = A\b;
err_rel = norm(sol1-sol2)/norm(sol1)

f = @() cauchySolver(x,y,R,S,b);
g = @() A\b;

timeit_cauchySolver = timeit(f)
timeit_backslash = timeit(g)

cputime_cauchySolver = cputime;
for k = 1:20
    sol1 = cauchySolver(x,y,R,S,b);
end
cputime_cauchySolver = (cputime - cputime_cauchySolver)/20

cputime_backslash = cputime;
for k = 1:20
    sol2 = A\b;
end
cputime_backslash = (cputime - cputime_backslash)/20

tic
for k = 1:20
    sol1 = cauchySolver(x,y,R,S,b);
end
elapsedtime_cauchySolver = toc/20

tic
for k = 1:20
    sol2 = A\b;
end
elapsedtime_backslash = toc/20