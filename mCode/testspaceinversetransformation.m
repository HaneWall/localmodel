%test_space for inverse fft 
central_wavelength = 800e-9;

c = 299792458;
amplitude = 1;
delay = 500e-15;
tau_int = 140e-15;
delta_t = 3e-18;
time = 0:delta_t:1000e-15;

omega = 2*pi * c / central_wavelength;
e_field = amplitude * exp(-2*log(2)*((time - delay)/tau_int).^2).*sin(omega.*(time - delay));
L = length(time);
n_fft = 2^nextpow2(L);

fft_e_field = fft(e_field, n_fft);
f = 1/delta_t*(0:(n_fft/2))/n_fft;
delta_f = 1/delta_t;
f_centra

