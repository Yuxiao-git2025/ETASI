function [eventOut, gridOut] = cal_ETASI_Evaluation(theta, time, mag0, smartSample)
% Evaluate fitted ETASI model at event times and on a plotting grid
[~, ~, details] = cal_ETASI_NegLogLik_Inner(theta, time, mag0, smartSample);

A     = theta(1);
c     = theta(2);
p     = theta(3);
alpha = theta(4);
mu    = theta(5);
b     = theta(6);
Tb    = theta(7);

N = length(time);
beta = b * log(10);

eventOut = struct();

eventOut.time = time;
eventOut.mag0 = mag0;

eventOut.R0 = details.R0_i;
eventOut.R = details.R_i;
eventOut.lambda_true = details.R0_i;
eventOut.lambda_observed = details.R_i;

eventOut.logRate = details.logRate;
eventOut.logMagDensity = details.logMagDensity;

eventOut.N0 = Tb * eventOut.R0;
eventOut.beta = beta;

eventOut.mu = mu * ones(N,1);

eventOut.trig = eventOut.R0 - mu;

eventOut.f0_GR = beta * exp(-beta * mag0);

eventOut.f_ETASI = exp(eventOut.logMagDensity);

eventOut.detection_rate_ratio = eventOut.R ./ eventOut.R0;

eventOut.Params = theta;

% Plotting grid
n_grid = 500;
T_grid = linspace(time(1), time(end), n_grid)';

R0_grid = zeros(n_grid,1);
R_grid = zeros(n_grid,1);
trig_grid = zeros(n_grid,1);

for k = 1:n_grid
    t = T_grid(k);
    idx = find(time < t);

    if isempty(idx)
        trig = 0;
    else
        dt = t - time(idx);
        trig = sum(A * exp(alpha * mag0(idx)) .* (c + dt).^(-p));
    end

    R0_grid(k) = mu + trig;
    R_grid(k) = cal_ETASI_Rate(R0_grid(k), Tb);
    trig_grid(k) = trig;
end

gridOut = struct();

gridOut.time = T_grid;
gridOut.R0 = R0_grid;
gridOut.R = R_grid;
gridOut.lambda_true = R0_grid;
gridOut.lambda_observed = R_grid;
gridOut.mu = mu * ones(n_grid,1);
gridOut.trig = trig_grid;
gridOut.N0 = Tb * R0_grid;
gridOut.detection_rate_ratio = R_grid ./ R0_grid;
gridOut.Params = theta;

end
