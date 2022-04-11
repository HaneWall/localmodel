function [value] = integrate_population_cb(ADK, delta_t, t)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    m = length(t);
    rho_sfi = zeros(size(t));
    %CN-Propagator
    for n = 1:(m-1)
        G = (1 - delta_t/2 * ADK(n))/(1 + delta_t/2 * ADK(n+1));
        E = delta_t/2 * (ADK(n+1) + ADK(n)) / (1 + delta_t/2 * ADK(n+1));
        rho_sfi(n+1) = G * rho_sfi(n) + E;
    end
    value = rho_sfi;
end

