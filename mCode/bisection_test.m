%bisection test

f = @(x) 3*x-3;

x_zero = bisection_err(f, -5, 5, 0.001);

function [z] = lin(x)
    z = 2*x-3;
end