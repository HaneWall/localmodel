function [rate] = ADK_rate(adk_0,amplitude, wavelength, time)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
rate = adk_0 * (sinusoidal_efield(amplitude, wavelength, time) / amplitude)^12;
end