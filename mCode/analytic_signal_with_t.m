function [z, phase] = analytic_signal_with_t(x, t, t_ref)
n_fft = length(x);
n_half = ceil(n_fft/2);
[~, idx_ref] = min(abs(t-t_ref)); 

% shift parameter --> get the reference time in the middle of the array
shift = n_half - idx_ref;
circ_shifted_x = circshift(x, shift);

% fft of circular shifted array 
X = fft(ifftshift(circ_shifted_x), n_fft);

% mask out the negative frequencies 
X_shifted = fftshift(X);
X_shifted(1:floor(n_fft/2)) = 0;
X_shifted(floor(n_fft/2)+1:end-1) = 2*X_shifted(floor(n_fft/2)+1:end-1);

% now we can inverse transform the masked signal 
z = ifft(X_shifted, n_fft);
z = ifftshift(z);

% reverse shift the envelope position
z = circshift(z, -shift);

% to get phase information we can work with X_shifted
X2 = X_shifted;
% mask out values that correspond to noise 
%X2(abs(X_shifted)<(max(X_shifted)/4)) = 0; 
phase = angle(X2);
end