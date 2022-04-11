% perp_kerr_vergleich.m
% x-probe, y-pump
clear all;

q = -1.60217662e-19;
me = 9.10938e-31;
n0 = 2.2e28; 
bandgaps = [7.5]; %in eV TODO: multiple bandgaps --> multiple adk0
wavelength_probe = 800e-9;
wavelength_pump = 2100e-9;
adk0 = 8.540703969149006e+12;

%varying pump intensities 
no_simulations = 30;
injection_first_harm = zeros(length(bandgaps),no_simulations);
brunel_first_harm = zeros(length(bandgaps),no_simulations);
kerr_first_harm = zeros(length(bandgaps),no_simulations);
overall_along_x_harm_perp = zeros(length(bandgaps),no_simulations);

% integration params
t_end = 1000e-15;
delta_t = 5e-18;
t = 0:delta_t:t_end;

%delay times
tau_pump = 430e-15;
tau_probe = 430e-15;

e_pump_ranges = linspace(3e16, 18e16, no_simulations);
amplitude_probe= intensity2amplitude(0.015e16); %0.015 TWcm^-2
amplitude_pump = intensity2amplitude(e_pump_ranges); %12 TWcm^-2

L = length(t);
n_fft = 2^nextpow2(L); %zero padding
f = 1/delta_t*(0:(n_fft/2))/n_fft;
delta_f = 1/delta_t;
f_pump = physconst('lightspeed') / 2.1e-06;
f_probe = physconst('lightspeed') / 8e-07;
first_harm = 2*f_pump + f_probe;
[~, idx] = min(abs(first_harm - f));
%idx = 893;

for b = 1:length(bandgaps)
    for i = 1:no_simulations
        amplitude_sum = amplitude_pump(i) + amplitude_probe;
        e_field_pump = gaussian_efield_new(amplitude_pump(i), wavelength_pump, 140e-15, tau_pump, t);
        e_field_probe = gaussian_efield_new(amplitude_probe, wavelength_probe, 45e-15, tau_probe, t);
        e_field = e_field_pump + e_field_probe;
        displacements_x_pump = displacement_x_new(bandgaps(b), max(e_field_pump), e_field_pump);
        displacements_x_probe = displacement_x_new(bandgaps(b), max(e_field_pump), e_field_probe);
        displacements_x_overall = displacement_x_new(bandgaps(b), max(e_field), e_field_probe);
        norm_e_field = sqrt(e_field_pump.^2 + e_field_probe.^2);
        ADK = ADK_rate_new(adk0, norm_e_field);
        rho_sfi = integrate_population_cb(ADK, delta_t, t);
        drho = cent_diff_n(rho_sfi, delta_t, 3);
        third_term_pump = cent_diff_n(displacements_x_pump.*drho, delta_t, 3);
        third_term_probe = cent_diff_n(displacements_x_overall.*drho, delta_t, 3);
        v0 = 0;
        brunel_current_density = n0 * q * q/me * e_field_probe.*rho_sfi;
        kerr_pump = kerr_current_density(e_field_pump, delta_t); 
        kerr_probe = kerr_current_density(e_field_probe, delta_t); 
        kerr = kerr_current_density(e_field, delta_t);
        %injection_current_density =  n0 * q * (third_term_pump + third_term_probe); 
        injection_current_probe = n0 * q * (third_term_probe); 
        injection_current_pump = n0 * q * (third_term_pump); 
        overall_current_x = brunel_current_density + 1/3*kerr + injection_current_probe;

        ft_brunel_current = fft(brunel_current_density, n_fft);
        ft_injection_current = fft(injection_current_probe, n_fft);
        ft_kerr_current = fft(kerr, n_fft);
        ft_overall_current_x = fft(overall_current_x, n_fft);
        
        P_brunel_current = abs(ft_brunel_current/n_fft).^2;
        P_injection_current = abs(ft_injection_current/n_fft).^2;
        P_kerr_current = abs(ft_kerr_current/n_fft).^2;
        P_overall_current_x = abs(ft_overall_current_x/n_fft).^2;
        if i == 20
            figure(1)
            semilogy(P_overall_current_x(1:n_fft/2 + 1),'black');
            hold on
            semilogy(P_kerr_current(1:n_fft/2 + 1), 'blue');
            semilogy(P_brunel_current(1:n_fft/2 + 1), 'red');
            semilogy(P_injection_current(1:n_fft/2 + 1), 'green');
            xline(idx);
        end

        brunel_first_harm(b, i) = P_brunel_current(idx);
        kerr_first_harm(b, i) = P_kerr_current(idx);
        injection_first_harm(b, i) = P_injection_current(idx);
        overall_along_x_harm_perp(b, i) = P_overall_current_x(idx);
    end
end

normed_brunel = brunel_first_harm./kerr_first_harm;
normed_kerr = kerr_first_harm./kerr_first_harm;
normed_injection = injection_first_harm./kerr_first_harm;

figure(2)
semilogy(e_pump_ranges, normed_kerr, 'black');
hold on
semilogy(e_pump_ranges, normed_brunel, 'b-.');
semilogy(e_pump_ranges, normed_injection, 'b-');

figure(3)
plot(e_pump_ranges, overall_along_x_harm_perp(1,:));