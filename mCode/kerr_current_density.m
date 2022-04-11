function [der_current_density] = kerr_current_density(e_field, delta_t)
    %UNTITLED7 Summary of this function goes here
    %   Detailed explanation goes here
    eps0 = 8.8541878128e-12;
    % third order chi SiO2
    chi_3 = 1.94e-22;
    p = eps0 * chi_3 * e_field.^3;
    current_dens = cent_diff_n(p, delta_t, 3);
    der_current_density = cent_diff_n(current_dens, delta_t, 3);
end



