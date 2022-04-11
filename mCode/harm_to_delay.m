%script that does compare the first harmonic n = 1 for 2*n omega_pump +
%omega_probe in respect to the delay (parallel polarization)
clear all;
%physical constants
q = -1.60217662e-19;
me = 9.10938e-31;
n0 = 2.2e28;                                     %molecular density for si02

%simulation parameters 
bandgap = -7.5; % in eV

amplitude_probe= intensity2amplitude(0.015e16);  %0.015 TWcm^-2
amplitude_pump = intensity2amplitude(12e16);     %12 TWcm^-2

wavelength_probe = 800e-9;
wavelength_pump = 2100e-9;

fwhm_pump = 140e-15;
fwhm_probe = 45e-15;
%adk0 = 25965409784120.4;
adk0 = 8.540703969149006e+12;

%integration parameters
t_end = 1000e-15;
delta_t = 5e-18;
t = 0:delta_t:t_end;
%n = 400000;
%t = linspace(0, t_end, n);
%delta_t = t_end/(n-1);

%delay times
delay_between_pulses = linspace(-80e-15, 80e-15, 40);
tau_pump = ones(1, 40) * 500e-15;
tau_probe = tau_pump + delay_between_pulses;

%allocate memory for first harmonic
L = length(t);
n_fft = 2^nextpow2(L); %zero padding
power_density_harmonic = zeros(length(delay_between_pulses), n_fft/2 + 1);

%getting the signal that is only produced by the pump pulse
% e_field_pump_solo = gaussian_efield_new(amplitude_pump, wavelength_pump, fwhm_pump, tau_pump(1), t);
% displacements_x_solo = displacement_x_new(7.5, max(e_field_pump_solo), e_field_pump_solo);
% ADK_solo = ADK_rate_new(adk0, e_field_pump_solo);
% rho_sfi_pump = integrate_population_cb(ADK_solo, delta_t, t);
% drho_solo = cent_diff_n(rho_sfi_pump, delta_t, 3);
% third_term_pump = cent_diff_n(displacements_x_solo.*drho_solo, delta_t, 3);
% v0 = 0;
% plasma_current_density_pump = n0 * q * (q/me * e_field_pump_solo.*rho_sfi_pump + v0*drho_solo+ third_term_pump);
% brunel_current_density_pump = n0 * q * q/me * e_field_pump_solo.*rho_sfi_pump;
% kerr_pump = kerr_current_density(e_field_pump_solo, delta_t); % unsure of prefactor 
% injection_current_density_pump =  n0 * q * third_term_pump; 
% overall_current_density_pump = plasma_current_density_pump + brunel_current_density_pump + kerr_pump + injection_current_density_pump;
% ft_overall_current_pump = fft(overall_current_density_pump, n_fft);
% whole_power_spec_pump = abs(ft_overall_current_pump/n_fft).^2;    
% power_density_harmonic_pump = whole_power_spec_pump(1:n_fft/2 + 1);

for i = 1:length(delay_between_pulses)
    e_field_pump = gaussian_efield_new(amplitude_pump, wavelength_pump, fwhm_pump, tau_pump(i), t);
    e_field_probe = gaussian_efield_new(amplitude_probe, wavelength_probe, fwhm_probe, tau_probe(i), t);
    e_field = e_field_pump + e_field_probe;
    displacements_x = displacement_x_new(bandgap, max(e_field), e_field);
    ADK = ADK_rate_new(adk0, e_field);
    
    rho_sfi = integrate_population_cb(ADK, delta_t, t);
    drho = cent_diff_n(rho_sfi, delta_t, 3);
    third_term = cent_diff_n(displacements_x.*drho, delta_t, 3);
    v0 = 0;
    plasma_current_density = n0 * q * (q/me * e_field.*rho_sfi + v0*drho + third_term);
    brunel_current_density = n0 * q * q/me * e_field.*rho_sfi;
    kerr = kerr_current_density(e_field, delta_t);  
    injection_current_density =  n0 * q * third_term; 
    overall_current_density = plasma_current_density + kerr;
    
    ft_plasma_current = fft(plasma_current_density,n_fft);
    ft_brunel_current = fft(brunel_current_density, n_fft);
    ft_injection_current = fft(injection_current_density, n_fft);
    ft_kerr_current = fft(kerr, n_fft); 
    ft_overall_current = fft(overall_current_density, n_fft);
    whole_power_spec = abs(ft_overall_current/n_fft).^2;    
    power_density_harmonic(i,:) = whole_power_spec(1:n_fft/2 + 1);
end
log_power_spec = log(power_density_harmonic);
imagesc(log_power_spec, [118 130]);
