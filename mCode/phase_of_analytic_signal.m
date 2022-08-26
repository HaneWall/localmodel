function [phase] = phase_of_analytic_signal(X)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
X2 = X;
threshold = max(abs(X)) / 2;
X2(abs(X)<1e-5) = 0;
phase = angle(X2);
end