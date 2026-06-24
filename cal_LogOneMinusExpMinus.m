function y = cal_LogOneMinusExpMinus(x)
% y = log(1 - exp(-x)), x > 0
if x <= 0
    y = -inf;
elseif x < 1e-6
    y = log(-expm1(-x));
else
    y = log1p(-exp(-x));
end
end
