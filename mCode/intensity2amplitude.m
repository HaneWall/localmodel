function [amplitude] = intensity2amplitude(intensity)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
c = 299792458;
eps0 = 8.8541878128e-12;
amplitude = sqrt(2/(c*eps0) * intensity);
end