function [der_current_density_vec] = kerr_current_binomial(efield_1, efield_2, delta_t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    eps0 = 8.8541878128e-12;
    % third order chi SiO2
    chi_3 = 1.94e-22;
    L = length(efield_1);
    p = zeros(3, L);
    p = eps0 .* chi_3 .* ((dot(efield_1, efield_1) + 2*dot(efield_1, efield_2) + dot(efield_2, efield_2)) .* (efield_1 + efield_2));
    current_dens_vec = zeros(3, L);
    der_current_density_vec = zeros(3, L);
    for l=1:3
        current_dens_vec(l,:) = gradient(p(l,:), delta_t);
        der_current_density_vec(l,:) = gradient(current_dens_vec(l,:) , delta_t);
    end
end

