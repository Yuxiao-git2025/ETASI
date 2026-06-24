function n=cal_ETASI_BranchingRatio(theta, T, magmax0)
% Branching ratio calculation
% theta = [A, c, p, alpha, mu, b, Tb]
% This is the branching ratio of the underlying classic ETAS process,
% before short-term incompleteness is applied
A     = theta(1);
c     = theta(2);
p     = theta(3);
alpha = theta(4);
b     = theta(6);
beta = b * log(10);
if abs(alpha - beta) < 1e-10 || abs(p - 1) < 1e-10
    n = NaN;
    return;
end
Naft = A/(1-p)/(alpha-beta) ...
    * (exp((alpha-beta)*magmax0) - 1) ...
    * ((c+T)^(1-p) - c^(1-p));

Nmag = -(1/beta) * (exp(-beta*magmax0) - 1);

n = Naft / Nmag;
end
