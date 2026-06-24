function [negLogL, badflag, details]=cal_ETASI_NegLogLik_Inner(theta, time, mag0, smartSample)
% Compute negative log-likelihood:
%   negLogL = int R(t)dt - sum log[f(mi,ti)] - sum log[R(ti)]
badflag = false;
details = struct();
A     = theta(1);
c     = theta(2);
p     = theta(3);
alpha = theta(4);
mu    = theta(5);
b     = theta(6);
Tb    = theta(7);

N = length(time);
beta = b * log(10);

if any(theta <= 0) || beta <= 0
    negLogL = 1e100;
    badflag = true;
    return;
end

R0_i = zeros(N,1);
R_i = zeros(N,1);
logRate = zeros(N,1);
logMagDensity = zeros(N,1);

for i = 1:N
    if i == 1
        trig = 0;
    else
        dt = time(i) - time(1:i-1);
        source = A * exp(alpha * mag0(1:i-1)) .* (c + dt).^(-p);
        trig = sum(source);
    end

    R0 = mu + trig;
    x = Tb * R0;

    if R0 <= 0 || x <= 0 || ~isfinite(R0)
        negLogL = 1e100;
        badflag = true;
        return;
    end

    R = cal_ETASI_Rate(R0, Tb);

    if R <= 0 || ~isfinite(R)
        negLogL = 1e100;
        badflag = true;
        return;
    end

    eMag = exp(-beta * mag0(i));
    logDen = cal_LogOneMinusExpMinus(x);

    logRate(i) = log(R);

    logMagDensity(i) = log(beta) + log(x) ...
        - beta*mag0(i) ...
        - x*eMag ...
        - logDen;

    R0_i(i) = R0;
    R_i(i) = R;
end
if any(~isfinite(logRate)) || any(~isfinite(logMagDensity))
    negLogL = 1e100;
    badflag = true;
    return;
end

R_sample = cal_RateSamples(theta, time, mag0, smartSample);
integral_R = trapz(smartSample, R_sample);

logL = sum(logMagDensity) + sum(logRate) - integral_R;
negLogL = -logL;

details.R0_i = R0_i;
details.R_i = R_i;
details.logRate = logRate;
details.logMagDensity = logMagDensity;
details.integral_R = integral_R;
details.logL = logL;

if ~isfinite(negLogL)
    negLogL = 1e100;
    badflag = true;
end

end
