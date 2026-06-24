function tau=cal_ETASI_TransformedTime(theta, time, mag0, smartSample)
% Time-rescaling transform:
%                   tau_i = int_0^{ti} R(u) du
% If the ETASI model is adequate, transformed events should follow a
% unit-rate Poisson process
R_sample = cal_RateSamples(theta, time, mag0, smartSample);
cumInt = zeros(length(smartSample),1);

for k = 2:length(smartSample)
    cumInt(k) = cumInt(k-1) + ...
        0.5 * (R_sample(k) + R_sample(k-1)) * ...
        (smartSample(k) - smartSample(k-1));
end
tau = interp1(smartSample, cumInt, time, 'linear', 'extrap');
tau = tau(:);
end
