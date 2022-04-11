function [rate] = ADK_rate_gaussian(adk_0,amplitude, wavelength, time, tau_int)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
rate = adk_0 * (gaussian_efield(amplitude, wavelength, tau_int, time) / amplitude)^12;
end