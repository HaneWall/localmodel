function [field] = gaussian_efield(amplitude, wavelength, tau_int, time)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
c = physconst('LightSpeed');
omega = 2*pi * c / wavelength;
field = amplitude * exp(-2*log(2)*((time - 3*tau_int)/tau_int).^2).*sin(omega.*time);
end