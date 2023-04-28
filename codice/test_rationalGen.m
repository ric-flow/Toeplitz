% n = 128;
% a = rand(n,1);
% b = rand(n,1);
% b(1) = a(1);
% T = toeplitz(a,b);
% % p = [3 -1 6];
% % q = [2 -5 1 -8];
% [p,q] = myPade_exp(6,6);
% [h,gamma,beta] = rationalConverter(p,q);
% A = polyvalm(q,T)\polyvalm(p,T);
% N = steinDisplace(A);
% f = @(c,r,y) toeplitz(c,r)\y;
% [G,B] = simpleRatGen(a,b,h,gamma,beta,f);
% % G = real(G);
% % B = real(B);
% norm(N - G*B')/norm(N)
%
% A = polyvalm(q,T)\polyvalm(p,T);
% N = steinDisplace(A);
% solver = @(c,r,poly,Y) polyvalm(poly,toeplitz(c,r))\Y;
% [G,B] = rationalGen(a,b,p,q,solver);
% norm(N - G*B')/norm(N)

L = 13;

err1 = zeros(1,L);
err2 = zeros(1,L);

[p,q] = myPade_exp(6,6);
theta = 0.55;

for k = 0:L-1
    n = 1024 + 256*k;
    [a,b,T] = Merton(n);
    a = a/1e6;
    b = b/1e6;
    T = T/1e6;

    f = @(a,b,rhs) toeplitz(a,b)\rhs;
    g = @(a,b,rhs) T2CSolver(a,b,rhs);

    [A1, dr1] = expmT(a,b,p,q,theta,f);
    [A2, dr2] = expmT(a,b,p,q,theta,g);
    A = expm(T);

    err1(k+1) = norm(A1-A,"fro")/norm(A,"fro");
    err2(k+1) = norm(A2-A,"fro")/norm(A,"fro");
end

err1 = log10(err1);
err2 = log10(err2);

figure(1)
plot(1024:256:4096,err1,'-*',1024:256:4096,err2,'-*')
title('Errore relativo')
xlabel('size')
ylabel('log_{10}(e.r.)')
legend({'backslash','GKO'})
