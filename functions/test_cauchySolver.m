% % test per il solver di sistemi cauchy
% % la nostra implementazione non è molto veloce
% % c'è comunque un guadagno nei tempi di cpu per matrici grandi (4000 o più)
% 
% n = 4000;
% dr = 4;
% b_size = 7;
% b = rand(n,b_size);
% R = rand(n,dr);
% S = rand(n,dr);
% x = 1:2:2*n;
% y = 2:2:2*n;
% A = cauchyConstructor(R,S,x,y);
% 
% sol1 = cauchySolver(x,y,R,S,b);
% sol2 = A\b;
% err_rel = norm(sol1-sol2)/norm(sol1)
% 
% f = @() cauchySolver(x,y,R,S,b);
% g = @() A\b;
% 
% timeit_cauchySolver = timeit(f)
% timeit_backslash = timeit(g)
% 
% cputime_cauchySolver = cputime;
% for k = 1:20
%     sol1 = cauchySolver(x,y,R,S,b);
% end
% cputime_cauchySolver = (cputime - cputime_cauchySolver)/20
% 
% cputime_backslash = cputime;
% for k = 1:20
%     sol2 = A\b;
% end
% cputime_backslash = (cputime - cputime_backslash)/20
% 
% tic
% for k = 1:20
%     sol1 = cauchySolver(x,y,R,S,b);
% end
% elapsedtime_cauchySolver = toc/20
% 
% tic
% for k = 1:20
%     sol2 = A\b;
% end
% elapsedtime_backslash = toc/20
% 
% clear

n = 8192;
size_rhs = 16;
% a = rand(n,1);
% b = rand(n,1);
% b(1) = a(1);
% T = toeplitz(a,b);
[a,b,T] = Merton(n);
rhs = rand(n,size_rhs);
x1 = T\rhs;
x2 = T2CSolver(a,b,rhs);

fT = fft([a;0;b(n:-1:2)]);
matvec = @(x) tMatVec2(fT,x);
M = sparse(diag(a(1)*ones(n,1))+diag(a(2)*ones(n-1,1),1)+diag(b(2)*ones(n-1,1),-1));

x3 =  gmresMultipleRhs(matvec,rhs,100,1e-12,4,M);

norm(x1-x2)/norm(x1)
norm(x1-x3)/norm(x1)

f = @() T\rhs;
g = @() T2CSolver(a,b,rhs);

h = @() gmresMultipleRhs(matvec,rhs,100,1e-12,4,M);

timeit(f)
timeit(g)
timeit(h)

% gmres con precondizionatore tridiagonale 
% è il più veloce per il Merton Model

function x = gmresnomsg(varargin)
    [x,~] = gmres(varargin{:});
end

function x = gmresMultipleRhs(mat,rhs,restart,tol,maxit,prec)
    x = zeros(size(rhs));
    for j = 1:size(rhs,2)
        x(:,j) = gmresnomsg(mat, rhs(:,j), restart, tol, maxit, prec);
    end
end