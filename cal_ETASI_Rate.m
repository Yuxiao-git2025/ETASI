function R=cal_ETASI_Rate(R0, Tb)
% R(t) = (1/Tb)*(1 - exp(-Tb*R0(t)))
x = Tb * R0;
R = -expm1(-x) ./ Tb;
end
