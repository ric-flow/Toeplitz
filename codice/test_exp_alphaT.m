% scalabilità

n = 2048;
[p,q] = myPade_exp(6,6);
theta = 0.55;

aux = (2.^(0:n-1))';
a = rand(n,1);
b = rand(n,1);
a(1) = -1;
b(1) = -1;
a = a./aux;
b = b./aux;

T = toeplitz(a,b);
figure(1)
plot(eig(T),'x')

t1 = zeros(1,9);
t2 = zeros(1,9);
t = zeros(1,9);
err1 = zeros(1,9);
err2 = zeros(1,9);

for k = 0:8
    if k > 0
        a = 2*a;
        b = 2*b;
    end
    solver1 = @(a,b,rhs) T2CSolver(a,b,rhs);
    solver2 = @(a,b,rhs) gmresStrang(a,b,rhs,100,1e-12,4);
    t1(k+1) = cputime;
    A1 = expmT(a,b,p,q,theta,solver1);
    t1(k+1) = cputime - t1(k+1);
    t2(k+1) = cputime;
    A2 = expmT(a,b,p,q,theta,solver2);
    t2(k+1) = cputime - t2(k+1);
    T = toeplitz(a,b);
    t(k+1) = cputime;
    A = expm(T);
    t(k+1) = cputime - t(k+1);
    err1(k+1) = norm(A-A1,"fro")/norm(A,"fro");
    err2(k+1) = norm(A-A2,"fro")/norm(A,"fro");
end

l = 0:8;
t = log2(t);
t1 = log2(t1);
t2 = log2(t2);

figure(2)
plot(l,t,'g-.*',l,t1,'b-.*',l,t2,'r-.*')
title('Confronto tempi di calcolo')
xlabel('log_2(alpha)')
ylabel('log_2(time)')
legend({'expm','expmt GKO','expmt gmres'})

err1 = log10(err1);
err2 = log10(err2);

figure(3)
plot(l,err1,'b-.*',l,err2,'r-.*')
title('Errore relativo')
xlabel('log_2(alpha)')
ylabel('log_{10}(e.r.)')
legend({'expmt GKO','expmt gmres'})


function x = gmresStrang(a,b,rhs,restart,tol,maxit)
    a = a(:);
    b = b(:);
    fT = fft([a;0;b(end:-1:2)]);
    matvec = @(x) tMatVec2(fT,x);
    c = strangPrec(a,b);
    c_aux = fft(c);
    precmatvec = @(x) ifft(fft(x)./c_aux);
    [x,~] = gmres(matvec,rhs,restart,tol,maxit,precmatvec);
end
