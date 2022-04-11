clear all;
%physical constants
q = -1.60217662e-19;
me = 0.635 * 9.10938e-31;
n0 = 2.2e28; %molecular density for si02
c = 299792458;

%simualtion parameters 
%amplitude_probe= intensity2amplitude(0.015e16); %0.015 TWcm^-2
amplitude_probe = 0; %0.015 TWcm^-2
amplitude_pump = intensity2amplitude(12e16); %12 TWcm^-2
amplitude_sum = amplitude_pump + amplitude_probe;
wavelength_probe = 800e-9;
wavelength_pump = 2100e-9;
%adk0 = 25965409784120.4;
adk0 = 8.540703969149006e+12;
t_end = 1000e-15;
delta_t = 5e-18;
t = 0:delta_t:t_end;

%delay times
tau_pump = 430e-15;
tau_probe = 430e-15;

%initiate arrays
e_field_pump = gaussian_efield_new(amplitude_pump, wavelength_pump, 140e-15, tau_pump, t);
e_field_probe = gaussian_efield_new(amplitude_probe, wavelength_probe, 45e-15, tau_probe, t);
e_field = e_field_pump + e_field_probe;
displacements_x = displacement_x_new(7.5, max(e_field), e_field);
ADK = ADK_rate_new(adk0, e_field);

rho_sfi = integrate_population_cb(ADK, delta_t, t);

%compute time derivative plasma current density 
drho = cent_diff_n(rho_sfi, delta_t, 3);
third_term = cent_diff_n(displacements_x.*drho, delta_t, 3);
v0 = 0;
plasma_current_density = n0 * q * (q/me * e_field.*rho_sfi + v0*drho+ third_term);
brunel_current_density = n0 * q * q/me * e_field.*rho_sfi;
kerr = kerr_current_density(e_field, delta_t); 
injection_current_density =  n0 * q * third_term; 
overall_current_density = plasma_current_density + kerr ;

figure()
%plot(t - tau_pump, e_field_pump.^2./max(e_field_pump.^2) .* max(rho_sfi),'.');
plot(t-tau_pump, adk0*(abs(e_field)./max(e_field)).^(13) ./max(adk0*(abs(e_field)./max(e_field)).^(13)))
hold on 
plot(t - tau_pump, e_field_probe.^2./max(e_field_probe.^2) .* max(rho_sfi), '-.');
hold on
plot(t - tau_pump, rho_sfi, '-');
figure()
plot(t, plasma_current_density);
figure()
L = length(t);
n_fft = 2^nextpow2(L); %zero padding
ft_plasma_current = fft(plasma_current_density,n_fft);
ft_brunel_current = fft(brunel_current_density, n_fft);
ft_injection_current = fft(injection_current_density, n_fft);
ft_kerr_current = fft(kerr, n_fft);
ft_overall_current = fft(overall_current_density, n_fft);

f = 1/delta_t*(0:(n_fft/2))/n_fft;
P_brunel_current = abs(ft_brunel_current/n_fft).^2;
P_injection_current = abs(ft_injection_current/n_fft).^2;
P_kerr_current = abs(ft_kerr_current/n_fft).^2;
P_overall_current = abs(ft_overall_current/n_fft).^2;

omegas = 2*pi.*f;
omega_pump = 2*pi * c / 2.1e-06;
omega_probe = 2*pi * c / 8e-07;
lambdas = 2*pi * c ./omegas;
semilogy(omegas,P_brunel_current(1:n_fft/2 + 1)/max(P_overall_current), 'red');
hold on
semilogy(omegas,P_injection_current(1:n_fft/2 + 1)/max(P_overall_current),'green');
hold on
semilogy(omegas,P_kerr_current(1:n_fft/2 + 1)/max(P_overall_current),'blue');
hold on
semilogy(omegas,P_overall_current(1:n_fft/2 + 1)/max(P_overall_current), 'black');
for k = 1:20
    %xline((2*k*omega_pump + omega_probe),'-');
    xline(k*omega_pump, '-.');
    %xline(k*omega_probe, ':');
end
%legend({'Brunel','Injection', 'Kerr', 'Overall', 'harmonics', 'pump', 'probe'})
legend({'Brunel','Injection', 'Kerr', 'Overall','pump'})
xlabel('\omega in rad/s')
ylabel('|FFT(\partial_t j)|^2')



