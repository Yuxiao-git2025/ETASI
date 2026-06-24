function negLogL = cal_ETASI_NegLoglik(z, time, mag0, LB, UB, smartSample)
% Negative log-likelihood of the ETASI model
theta = cal_UnitToParams(z, LB, UB);
[negLogL, badflag] = cal_ETASI_NegLogLik_Inner(theta, time, mag0, smartSample);
if badflag || ~isfinite(negLogL)
    negLogL = 1e20;
end

end
