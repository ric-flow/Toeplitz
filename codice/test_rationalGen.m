n = 1024;
a = rand(n,1);
b = rand(n,1);
b(1) = a(1);
T = toeplitz(a,b);
% p = [3 -1 6];
% q = [2 -5 1 -8];
[p,q] = myPade_exp(6,6);
[h,gamma,beta] = rationalConverter(p,q);
A = rationalEval(T,h,gamma,beta);
N = steinDisplace(A);
f = @(c,r,y) toeplitz(c,r)\y;
[G,B] = simpleRatGen(a,b,h,gamma,beta,f);
% G = real(G);
% B = real(B);
norm(N - G*B')/norm(N)