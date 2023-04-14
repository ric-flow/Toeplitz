function [p,q] = myPade_exp(k,s)

% genera l'approssimante di Padé di grado (k,s) della funzione esponenziale
% INPUT:
% k : grado del polinomio numeratore
% s : grado del polinomio denominatore
% OUTPUT:
% p : polinomio numeratore
% q : polinomio denominatore
% NOTA: richiede MATLAB symbolic toolbox

syms z;
syms psym;
syms qsym;
f = pade(exp(z),'Order',[k,s]);

if k > s
    psym = children(f,2);
    c = children(f,3);
    psym = c*psym;
    qsym = children(f,1);
    qsym = 1/qsym;
end

if k < s
    psym = children(f,1);
    c = children(f,3);
    psym = c*psym;
    qsym = children(f,2);
    qsym = 1/qsym;
end

if k == s && mod(k,2) == 1
    psym = children(f,2);
    c = children(f,3);
    psym = c*psym;
    qsym = children(f,1);
    qsym = 1/qsym;
end

if k == s && mod(k,2) == 0
    psym = children(f,2);
    qsym = children(f,1);
    qsym = 1/qsym;
end

p = sym2poly(psym);
q = sym2poly(qsym);

end