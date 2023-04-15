% test expmT

n = 2048;

[a,b,T] = Merton(n);

f = @(a,b,rhs) toeplitz(a,b)\rhs;
g = @(a,b,rhs) T2CSolver(a,b,rhs);
h = @(a,b,rhs) gmresMerton(a,b,rhs,100,1e-11,4);

[p,q] = myPade_exp(6,6);
theta = 0.55;

A = expm(T);
A1 = expmT(a,b,p,q,theta,f);
A2 = expmT(a,b,p,q,theta,g);
A3 = expmT(a,b,p,q,theta,h);

err1 = norm(A-A1,"fro")/norm(A,"fro")
err2 = norm(A-A2,"fro")/norm(A,"fro")
err3 = norm(A-A3,"fro")/norm(A,"fro")

esp = @() expm(T);
exp1 = @() expmT(a,b,p,q,theta,f);
exp2 = @() expmT(a,b,p,q,theta,g);
exp3 = @() expmT(a,b,p,q,theta,h);

timeit(esp)
timeit(exp1)
timeit(exp2)
timeit(exp3)

% gmres con precondizionatore tridiagonale 
% è il più veloce per il Merton Model

function x = gmresMerton(a,b,rhs,restart,tol,maxit)
    a = a(:);
    b = b(:);
    n = length(a);
    fT = fft([a;0;b(n:-1:2)]);
    matvec = @(x) tMatVec2(fT,x);
    prec = sparse(diag(a(1)*ones(n,1))+diag(a(2)*ones(n-1,1),1)+diag(b(2)*ones(n-1,1),-1));
    [x,~] = gmres(matvec,rhs,restart,tol,maxit,prec);
end
