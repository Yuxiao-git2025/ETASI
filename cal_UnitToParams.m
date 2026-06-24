function theta=cal_UnitToParams(z, LB, UB)
% Convert normalized variables z in [0,1] to true ETASI parameters
% theta = [A, c, p, alpha, mu, b, Tb]
z = z(:)';
A     = 10^(UB(1)*z(1) + LB(1)*(1-z(1)));
c     = 10^(UB(2)*z(2) + LB(2)*(1-z(2)));
p     = UB(3)*z(3) + LB(3)*(1-z(3));
alpha = UB(4)*z(4) + LB(4)*(1-z(4));
mu    = 10^(UB(5)*z(5) + LB(5)*(1-z(5)));
b     = UB(6)*z(6) + LB(6)*(1-z(6));
Tb    = 10^(UB(7)*z(7) + LB(7)*(1-z(7)));

theta = [A, c, p, alpha, mu, b, Tb];
end
