function [e_field] = gaussian_efield_new(amplitude, central_wavelength, tau_int, delay, time)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
c = 299792458;
omega = 2*pi * c / central_wavelength;
e_field = amplitude * exp(-2*log(2)*((time - delay)/tau_int).^2).*sin(omega.*time);
end