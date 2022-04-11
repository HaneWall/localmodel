function [x_0] = displacement_x_new(bandgap_in_ev, normed_field, e_field)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
q = 1.60217662e-19;
bandgap_in_j = bandgap_in_ev * q;
%x_0 = e_field .* bandgap_in_j ./(q* max_amplitude.^2);
x_0 = e_field .* bandgap_in_j ./(q* normed_field.^2);
end