% Modeling:
%   R0(t) = mu + sum_{tj<t} A*exp(alpha*(mj-Mmin))*(c+t-tj)^(-p)
% 
% Apparent ETASI rate:
%   R(t) = (1/Tb) * (1 - exp(-Tb*R0(t)))
%
% Apparent magnitude density:
%   f(m,t) = beta*N0(t)*exp(-beta*(m-Mmin)) ...
%            * exp(-N0(t)*exp(-beta*(m-Mmin))) ...
%            / (1 - exp(-N0(t)))
% where:
%   beta  = b*log(10)
%   N0(t) = Tb*R0(t)
%
% Inputs:
%   time   : event times, preferably in days
%   mag    : event magnitudes
%   Mmin   : completeness magnitude Mc
%   theta0 : initial ETASI parameters
%            [A, c, p, alpha, mu, b, Tb]
%
% Output:
%   Result : structure containing fitted parameters, likelihood, AIC,
%            event-wise rates, transformed time, and model settings
function Result=cal_ETASI(time, mag, Mmin, Param0)
% -------------------------------------------------------------------------
% Basic setting
% -------------------------------------------------------------------------
time = time(:);
mag = mag(:);
[time, sortId] = sort(time);
mag = mag(sortId);
N = length(time);

if length(Param0) ~= 7
    error('Parameters must be [A, c, p, alpha, mu, b, Tb]');
end

% ETASI fitting uses relative magnitudes
mag0=mag-Mmin;

% Shift time to start from zero, following the Python code.
time0 = time - time(1);
T = time0(end);
magmax0 = max(mag0);

% -------------------------------------------------------------------------
% Parameter bounds
% -------------------------------------------------------------------------
% Normalized optimization variable z is in [0,1]
%
% Log-scale parameters:
%   A, c, mu, Tb
%
% Linear-scale parameters:
%   p, alpha, b

min_A  = log10(1e-6);
max_A  = log10(5e1);

min_c  = log10(1e-6);
max_c  = log10(1e0);

min_p  = 0.2;
max_p  = 3.0;

min_al = 0.1;
max_al = 5.5;

min_mu = log10(1e-4);
max_mu = log10(100);

min_b  = 0.2;
max_b  = 2.0;

min_Tb = log10(1e-6);
max_Tb = log10(1/24);

LB = [min_A, min_c, min_p, min_al, min_mu, min_b, min_Tb];
UB = [max_A, max_c, max_p, max_al, max_mu, max_b, max_Tb];
Param0 = Param0(:)';

z0 = cal_ParamsToUnit(Param0, LB, UB);
z0 = min(max(z0, 0), 1);

% -------------------------------------------------------------------------
% Smart samples for numerical integration
% -------------------------------------------------------------------------
smartSample = cal_BuildSamples(time0);

% -------------------------------------------------------------------------
% Initial diagnostics
% -------------------------------------------------------------------------
n0 = cal_ETASI_BranchingRatio(Param0, T, magmax0);
fprintf('\n=== ETASI INITIAL PARAMETERS ===\n');
fprintf('A      = %.10g\n', Param0(1));
fprintf('c      = %.10g\n', Param0(2));
fprintf('p      = %.10g\n', Param0(3));
fprintf('alpha  = %.10g\n', Param0(4));
fprintf('mu     = %.10g\n', Param0(5));
fprintf('b      = %.10g\n', Param0(6));
fprintf('Tb     = %.10g\n', Param0(7));
fprintf('n      = %.10g\n', n0);

% -------------------------------------------------------------------------
% Optimization
% -------------------------------------------------------------------------
objfun = @(z) cal_ETASI_NegLoglik(z, time0, mag0, LB, UB, smartSample);
if exist('fmincon', 'file')==2
    options = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'Algorithm', 'interior-point', ...
        'MaxIterations', 300, ...
        'MaxFunctionEvaluations', 3e5, ...
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-6);

    [z_opt, negLogL, exitflag, output] = fmincon( ...
        objfun, z0, [], [], [], [], zeros(7,1), ones(7,1), [], options);
else
    options = optimset( ...
        'Display', 'iter', ...
        'MaxIter', 300, ...
        'MaxFunEvals', 3e5, ...
        'TolX', 1e-6, ...
        'TolFun', 1e-6);

    boundedObj = @(z) cal_ETASI_Bounds(z, time0, mag0, LB, UB, smartSample);
    [z_opt, negLogL, exitflag, output] = fminsearch(boundedObj, z0, options);
    z_opt = min(max(z_opt, 0), 1);
    negLogL = cal_ETASI_NegLoglik(z_opt, time0, mag0, LB, UB, smartSample);
end

theta_opt = cal_UnitToParams(z_opt, LB, UB);

A_opt     = theta_opt(1);
c_opt     = theta_opt(2);
p_opt     = theta_opt(3);
alpha_opt = theta_opt(4);
mu_opt    = theta_opt(5);
b_opt     = theta_opt(6);
Tb_opt    = theta_opt(7);

logL = -negLogL;
Lenparams = 7;
AIC = 2*Lenparams - 2*logL;

n_opt = cal_ETASI_BranchingRatio(theta_opt, T, magmax0);

% -------------------------------------------------------------------------
% Evaluate fitted model
% -------------------------------------------------------------------------
[eventOut, gridOut] = cal_ETASI_Evaluation(theta_opt, time0, mag0, smartSample);
TransN = cal_ETASI_TransformedTime(theta_opt, time0, mag0, smartSample);
% -------------------------------------------------------------------------
% Save results
% -------------------------------------------------------------------------
Result = struct();

Result.time = time;
Result.time0 = time0;
Result.mag = mag;
Result.mag0 = mag0;
Result.Mmin = Mmin;

Result.theta0 = Param0;
Result.theta = theta_opt;

Result.A = A_opt;
Result.c = c_opt;
Result.p = p_opt;
Result.alpha = alpha_opt;
Result.mu = mu_opt;
Result.b = b_opt;
Result.Tb = Tb_opt;

Result.Params_name = {'A', 'c', 'p', 'alpha', 'mu', 'b', 'Tb'};
Result.Params_final = theta_opt;

Result.negLogL = negLogL;
Result.logL = logL;
Result.AIC = AIC;
Result.exitflag = exitflag;
Result.output = output;

Result.branching_ratio = n_opt;

Result.bounds_LB = LB;
Result.bounds_UB = UB;
Result.z_opt = z_opt(:);

Result.event = eventOut;
Result.grid = gridOut;

Result.TransN = TransN;
Result.NumEvents = (1:N)';

fprintf('\n=== ETASI FIT FINISHED ===\n');
fprintf('A      = %.6g\n', A_opt);
fprintf('c      = %.6g\n', c_opt);
fprintf('p      = %.6g\n', p_opt);
fprintf('alpha  = %.6g\n', alpha_opt);
fprintf('mu     = %.6g\n', mu_opt);
fprintf('b      = %.4g\n', b_opt);
fprintf('Tb     = %.4g\n', Tb_opt);
fprintf('n      = %.6g\n', n_opt);
fprintf('logL   = %.6g\n', logL);
fprintf('AIC    = %.6g\n', AIC);
fprintf('Exitflag = %d\n', exitflag);

end
