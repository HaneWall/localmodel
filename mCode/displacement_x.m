function [x_0] = displacement_x(bandgap_in_ev, amplitude, time, wavelength)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
q = 1.60217662*10^(-19);
bandgap_in_j = bandgap_in_ev * q;
x_0 = sinusoidal_efield(amplitude, wavelength, time) * bandgap_in_j /(q * amplitude^2);
end