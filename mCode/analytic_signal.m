function [z, phase] = analytic_signal(x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n_fft = length(x);
X = fft(ifftshift(x), n_fft);

% mask out the negative frequencies 
X_shifted = fftshift(X);
X_shifted(1:floor(n_fft/2)) = 0;
X_shifted(floor(n_fft/2)+1:end-1) = 2*X_shifted(floor(n_fft/2)+1:end-1);

% now we can inverse transform the masked signal 
z = ifft(X_shifted, n_fft);
z = ifftshift(z);

% to get phase information we can work with X_shifted
X2 = X_shifted;
% mask out values that correspond to noise
X2(abs(X_shifted)<(max(X_shifted)/4)) = 0; 
phase = angle(X2);
end