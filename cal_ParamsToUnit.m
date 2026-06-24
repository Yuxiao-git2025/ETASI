function z = cal_ParamsToUnit(theta, LB, UB)
% Convert true ETASI parameters to normalized variables in [0,1].
%
% theta = [A, c, p, alpha, mu, b, Tb]

A     = theta(1);
c     = theta(2);
p     = theta(3);
alpha = theta(4);
mu    = theta(5);
b     = theta(6);
Tb    = theta(7);

z = zeros(1,7);

z(1) = (log10(A)  - LB(1)) / (UB(1) - LB(1));
z(2) = (log10(c)  - LB(2)) / (UB(2) - LB(2));
z(3) = (p         - LB(3)) / (UB(3) - LB(3));
z(4) = (alpha     - LB(4)) / (UB(4) - LB(4));
z(5) = (log10(mu) - LB(5)) / (UB(5) - LB(5));
z(6) = (b         - LB(6)) / (UB(6) - LB(6));
z(7) = (log10(Tb) - LB(7)) / (UB(7) - LB(7));

z = z(:);
end
