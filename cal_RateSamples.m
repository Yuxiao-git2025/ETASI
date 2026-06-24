function R_sample = cal_RateSamples(theta, time, mag0, sampleTime)
% Evaluate apparent ETASI rate R(t) on arbitrary sample times
A     = theta(1);
c     = theta(2);
p     = theta(3);
alpha = theta(4);
mu    = theta(5);
Tb    = theta(7);

R_sample = zeros(length(sampleTime),1);

for k = 1:length(sampleTime)
    t = sampleTime(k);
    idx = find(time < t);

    if isempty(idx)
        R0 = mu;
    else
        dt = t - time(idx);
        trig = sum(A * exp(alpha * mag0(idx)) .* (c + dt).^(-p));
        R0 = mu + trig;
    end

    R_sample(k) = cal_ETASI_Rate(R0, Tb);
end

end
