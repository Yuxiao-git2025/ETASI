function obj=cal_ETASI_Bounds(z, time, mag0, LB, UB, smartSample)
% Bounded objective for fminsearch fallback
if any(z < 0) || any(z > 1) || any(~isfinite(z))
    obj = 1e20 + 1e10*sum((min(z,0)).^2 + (max(z-1,0)).^2);
    return;
end
obj = Cost_ETASI(z, time, mag0, LB, UB, smartSample);
end
