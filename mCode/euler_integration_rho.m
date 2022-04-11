clear all;
amplitude = intensity2amplitude(5e16); %14 TWcm^2
n0 = 2.3e28;
q = -1.60217662*10^(-19);
me = 9.10938*10^(-31);


wavelength = 800*10^(-9);
adk0 = 2.596429371763225e+13;
t_end = 60*10^(-15);
rho0 = 0;
n = 5000;
t = linspace(0, t_end, n);
delta_t = t_end/(n-1);
e_field = gaussian_efield(amplitude, wavelength, 10e-15, t);
displacements_x = displacement_x_new(7.5, amplitude, e_field);
ADK = ADK_rate_new(adk0, e_field);

rho_sfi_euler = integrate_population_cb(ADK, delta_t, t);

%compute time derivative plasma current density 
drho = cent_diff_n(rho_sfi_euler, delta_t, 3);
third_term = cent_diff_n(displacements_x.*drho, delta_t, 3);
v0 = 0;

brunel_current_density = n0 * q *(q/me * e_field.*rho_sfi_euler);
kerr_current_density = kerr_current_density(e_field, delta_t); % unsure of prefactor 
ionization_current_density = n0 * q * third_term;
plasma_current_density = brunel_current_density + kerr_current_density + ionization_current_density;
figure()
plot(t, rho_sfi_euler, '-')
hold on
plot(t, abs(e_field)./max(e_field));
figure()
plot(t, plasma_current_density)
figure()
L = length(t);
n_fft = 2^nextpow2(L); %zero_padding
ff_plasma_current_density = fft(plasma_current_density, n_fft);
ft_brunel_current_density = fft(brunel_current_density, n_fft);
ft_ionization_current_density = fft(ionization_current_density, n_fft);
ft_kerr_current = fft(kerr_current_density, n_fft);
f = 1/delta_t*(0:(n_fft/2))/n_fft;
P_plasma_current = abs(ff_plasma_current_density/n_fft).^2;
P_brunel_current = abs(ft_brunel_current_density/n_fft).^2;
P_ionization_current = abs(ft_ionization_current_density/n_fft).^2;
P_kerr_current = abs(ft_kerr_current/n_fft).^2;

omegas = 2*pi.*f;
omega_0 = 2*pi*physconst('lightspeed')/wavelength;
semilogy(omegas/omega_0,P_plasma_current(1:n_fft/2 + 1)./max(P_plasma_current),'black') 
hold on 
semilogy(omegas/omega_0,P_ionization_current(1:n_fft/2 + 1)./max(P_plasma_current),'red')
semilogy(omegas/omega_0,P_brunel_current(1:n_fft/2 + 1)./max(P_plasma_current), 'green')
semilogy(omegas/omega_0,P_kerr_current(1:n_fft/2 + 1)./max(P_plasma_current), 'blue')
legend({'overall', 'injection', 'brunel', 'kerr'});
for k = 1:20
    if mod(k, 2) == 0
        xline(k, '-.')
    else
        xline(k)
    end
end
xlabel('Harmonische')
ylabel('|FFT(\partial_t j)|^2')



