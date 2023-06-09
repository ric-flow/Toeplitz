% % test expmT
% 
% t = zeros(1,12);
% t1 = zeros(1,12);
% t2 = zeros(1,12);
% t3 = zeros(1,12);
% err1 = zeros(1,12);
% err2 = zeros(1,12);
% err3 = zeros(1,12);
% 
% for k = 1:12
% 
% n = 512*k;
% 
% [a,b,T] = Merton(n);
% 
% f = @(a,b,rhs) toeplitz(a,b)\rhs;
% g = @(a,b,rhs) T2CSolver(a,b,rhs);
% h = @(a,b,rhs) gmresMerton(a,b,rhs,100,1e-11,4);
% 
% [p,q] = myPade_exp(6,6);
% theta = 0.55;
% 
% t(k) = cputime;
% A = expm(T);
% t(k) = cputime - t(k);
% t1(k) = cputime;
% A1 = expmT(a,b,p,q,theta,f);
% t1(k) = cputime - t1(k);
% t2(k) = cputime;
% A2 = expmT(a,b,p,q,theta,g);
% t2(k) = cputime - t2(k);
% t3(k) = cputime;
% A3 = expmT(a,b,p,q,theta,h);
% t3(k) = cputime - t3(k);
% 
% err1(k) = norm(A-A1,"fro")/norm(A,"fro");
% err2(k) = norm(A-A2,"fro")/norm(A,"fro");
% err3(k) = norm(A-A3,"fro")/norm(A,"fro");
% 
% % esp = @() expm(T);
% % exp1 = @() expmT(a,b,p,q,theta,f);
% % exp2 = @() expmT(a,b,p,q,theta,g);
% % exp3 = @() expmT(a,b,p,q,theta,h);
% % 
% % timeit(esp)
% % timeit(exp1)
% % timeit(exp2)
% % timeit(exp3)
% 
% end
% 
% figure(1)
% plot(1:12,t,1:12,t1,1:12,t2,1:12,t3)
% 
% figure(2)
% plot(1:12,err1,1:12,err2,1:12,err3)
% 
% % gmres con precondizionatore tridiagonale 
% % è il più veloce per il Merton Model



n = 2^14;
a_aux = rand(n,1);
a_aux(1) = - 1;
b_aux = rand(n,1);
b_aux(1) = a_aux(1);
a_aux = 16 * a_aux;
b_aux = 16 * b_aux;

[p,q] = myPade_exp(6,6);
theta = 0.55;


% test matrice con decadimento polinomiale

aux = ((1:n).^4)';
a_aux1 = a_aux./aux;
b_aux1 = b_aux./aux;
T = toeplitz(a_aux1(1:2000),b_aux1(1:2000));

figure(1)
plot(eig(T),'x');

t1 = zeros(1,7);
t2 = zeros(1,7);
t = zeros(1,6);
err1 = zeros(1,6);
err2 = zeros(1,6);

for k = 8:14
    a = a_aux1(1:2^k);
    b = b_aux1(1:2^k);
    solver1 = @(a,b,rhs) T2CSolver(a,b,rhs);
    solver2 = @(a,b,rhs) gmresStrang(a,b,rhs,100,1e-12,4);
    t1(k-7) = cputime;
    A1 = expmT(a,b,p,q,theta,solver1);
    t1(k-7) = cputime - t1(k-7);
    t2(k-7) = cputime;
    A2 = expmT(a,b,p,q,theta,solver2);
    t2(k-7) = cputime - t2(k-7);
    if k < 14
        T = toeplitz(a,b);
        t(k-7) = cputime;
        A = expm(T);
        t(k-7) = cputime - t(k-7);
        err1(k-7) = norm(A-A1,"fro")/norm(A,"fro");
        err2(k-7) = norm(A-A2,"fro")/norm(A,"fro");
    end
end
l = 9:14;
t = log2(t(2:end));
t1 = log2(t1(2:end));
t2 = max(log2(t2(2:end)),-8);

figure(2)
plot(9:13,t,'-.*',l,t1,'-.*',l,t2,'-.*')
title('Confronto tempi di calcolo')
xlabel('log_2(size)')
ylabel('log_2(time)')
legend({'expm','GKO','prec. gmres'})

err1 = log10(err1(2:end));
err2 = log10(err2(2:end));

figure(3)
plot(9:13,err1,'-*',9:13,err2,'-*')
title('Errore relativo')
xlabel('log_2(size)')
ylabel('log_{10}(e.r.)')
legend({'GKO','prec. gmres'})


[~,dr] = expmT(a,b,p,q,theta,solver2);
figure(4)
plot(0:length(dr)-1,dr,'-o')
title('Crescita displacement rank nella fase di squaring')
xlabel('passo di squaring')
ylabel('displacement rank')


% test matrice con decadimento esponenziale

aux = (2.^(1:n))';
a_aux2 = a_aux./aux;
b_aux2 = b_aux./aux;

T = toeplitz(a_aux2(1:2000),b_aux2(1:2000));
figure(5)
plot(eig(T),'x')

t1 = zeros(1,7);
t2 = zeros(1,7);
t = zeros(1,6);
err1 = zeros(1,6);
err2 = zeros(1,6);

for k = 8:14
    a = a_aux2(1:2^k);
    b = b_aux2(1:2^k);
    solver1 = @(a,b,rhs) T2CSolver(a,b,rhs);
    solver2 = @(a,b,rhs) gmresStrang(a,b,rhs,100,1e-12,4);
    t1(k-7) = cputime;
    A1 = expmT(a,b,p,q,theta,solver1);
    t1(k-7) = cputime - t1(k-7);
    t2(k-7) = cputime;
    A2 = expmT(a,b,p,q,theta,solver2);
    t2(k-7) = cputime - t2(k-7);
    if k < 14
        T = toeplitz(a,b);
        t(k-7) = cputime;
        A = expm(T);
        t(k-7) = cputime - t(k-7);
        err1(k-7) = norm(A-A1,"fro")/norm(A,"fro");
        err2(k-7) = norm(A-A2,"fro")/norm(A,"fro");
    end
end
l = 9:14;
t = log2(t(2:end));
t1 = log2(t1(2:end));
t2 = max(log2(t2(2:end)),-8);

figure(6)
plot(9:13,t,'-.*',l,t1,'-.*',l,t2,'-.*')
title('Confronto tempi di calcolo')
xlabel('log_2(size)')
ylabel('log_2(time)')
legend({'expm','GKO','prec. gmres'})

err1 = log10(err1(2:end));
err2 = log10(err2(2:end));

figure(7)
plot(9:13,err1,'-*',9:13,err2,'-*')
title('Errore relativo')
xlabel('log_2(size)')
ylabel('log_{10}(e.r.)')
legend({'GKO','prec. gmres'})

[~,dr] = expmT(a,b,p,q,theta,solver2);
figure(8)
plot(0:length(dr)-1,dr,'-o')
title('Crescita displacement rank nella fase di squaring')
xlabel('passo di squaring')
ylabel('displacement rank')



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


% function x = gmresMerton(a,b,rhs,restart,tol,maxit)
%     a = a(:);
%     b = b(:);
%     n = length(a);
%     fT = fft([a;0;b(n:-1:2)]);
%     matvec = @(x) tMatVec2(fT,x);
%     prec = sparse(diag(a(1)*ones(n,1))+diag(a(2)*ones(n-1,1),1)+diag(b(2)*ones(n-1,1),-1));
%     [x,~] = gmres(matvec,rhs,restart,tol,maxit,prec);
% end
