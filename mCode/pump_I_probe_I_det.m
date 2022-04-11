% script parallel pump/probe polarization
% E1 = [E_probe, 0, 0], E2 = [0, E_pump, 0], Detector=[0, 1, 0]
clear all;
q = -1.60217662e-19;
me = 9.10938e-31;
n0 = 2.2e28; 
bandgaps = [7.5]; %in eV TOO: multiple bandgaps --> multiple adk0
wavelength_probe = 800e-9;
wavelength_pump = 2100e-9;
adk0 = 8.540703969149006e+12;

% integration params
t_end = 1000e-15;
delta_t = 5e-18;
t = 0:delta_t:t_end;
L = length(t);

%varying pump intensities 
no_simulations = 30;
injection_first_harm = zeros(length(bandgaps),no_simulations);
brunel_first_harm = zeros(length(bandgaps),no_simulations);
kerr_first_harm = zeros(length(bandgaps),no_simulations);
overall_along_y_harm = zeros(length(bandgaps),no_simulations);

%delay times
tau_pump = 430e-15;
tau_probe = 430e-15;

e_pump_ranges = linspace(3e16, 18e16, no_simulations);
amplitude_probe= intensity2amplitude(0.015e16); %0.015 TWcm^-2
amplitude_pump = intensity2amplitude(e_pump_ranges);
e_field_probe = zeros(3, L);
e_field_pump = zeros(3, L);
third_term = zeros(3, L);
e_field_probe(1,:) = gaussian_efield_new(amplitude_probe, wavelength_probe, 45e-15, tau_probe, t);


n_fft = 2^nextpow2(L); %zero padding
f = 1/delta_t*(0:(n_fft/2))/n_fft;
delta_f = 1/delta_t;
f_pump = physconst('lightspeed') / 2.1e-06;
f_probe = physconst('lightspeed') / 8e-07;
first_harm = 2*f_pump + f_probe;
[~, idx] = min(abs(first_harm - f));
[~, idx_pump] = min(abs(f_pump - f));
[~, idx_probe] = min(abs(f_probe - f));

for b = 1:length(bandgaps)
    for i = 1:no_simulations
        amplitude_sum = amplitude_pump(i) + amplitude_probe;
        e_field_pump(2,:) = gaussian_efield_new(amplitude_pump(i), wavelength_pump, 140e-15, tau_pump, t);
        e_field = e_field_pump + e_field_probe;
        displacements_x = displacement_x_new(bandgaps(b), max(max(e_field)), e_field);
        for j = 1:L
            normed_e_field(:,j) = norm(e_field(:,j));
        end
        ADK = ADK_rate_new(adk0, normed_e_field);
        rho_sfi = integrate_population_cb(ADK, delta_t, t);
        drho = cent_diff_n(rho_sfi, delta_t, 3);
        third_term(2,:) = cent_diff_n(displacements_x(2,:).*drho, delta_t, 3);
        v0 = 0;
        brunel_current_density = n0 * q * q/me * e_field.*rho_sfi;
        kerr = kerr_current_density_vec(e_field_pump, e_field_probe, delta_t); 
        injection_current_density =  n0 .* q .* third_term; 
        overall_current_y = brunel_current_density(2,:) + kerr(2,:) + injection_current_density(2,:);

        ft_brunel_current = fft(brunel_current_density(2,:), n_fft);
        ft_injection_current = fft(injection_current_density(2,:), n_fft);
        ft_kerr_current = fft(kerr(2,:), n_fft);
        ft_overall_current_y = fft(overall_current_y, n_fft);
        P_brunel_current = abs(ft_brunel_current/n_fft).^2;
        P_injection_current = abs(ft_injection_current/n_fft).^2;
        P_kerr_current = abs(ft_kerr_current/n_fft).^2;
        P_overall_current = abs(ft_overall_current_y/n_fft).^2;
        brunel_first_harm(b, i) = P_brunel_current(idx);
        kerr_first_harm(b, i) = P_kerr_current(idx);
        injection_first_harm(b, i) = P_injection_current(idx);
        overall_along_y_harm(b, i) = P_overall_current(idx);
        if i == 20
            figure(1)
            semilogy(2*pi*f, P_overall_current(1:n_fft/2 + 1),'black');
            hold on
            semilogy(2*pi*f, P_kerr_current(1:n_fft/2 + 1), 'blue');
            semilogy(2*pi*f, P_brunel_current(1:n_fft/2 + 1), 'red');
            semilogy(2*pi*f, P_injection_current(1:n_fft/2 + 1), 'green');
            xline(2*pi*delta_f/n_fft*idx);
            xline(2*pi*delta_f/n_fft*idx_pump, '-.');
            xline(2*pi*delta_f/n_fft*idx_probe, '--');
            xlabel('$\omega$ in $rads^{-1}$','interpreter','latex');
            legend('Overall', 'Kerr', 'Brunel', 'Injection', '$2\omega_{pump} + \omega{_probe}$', '$\omega_{pump}$', '$\omega_{probe}$', 'Interpreter','latex');
        end
    end
end

normed_brunel = brunel_first_harm./kerr_first_harm;
normed_kerr = kerr_first_harm./kerr_first_harm;
normed_injection = injection_first_harm./kerr_first_harm;

figure(2)
plot(e_pump_ranges, overall_along_y_harm(1,:))


figure(3)
semilogy(e_pump_ranges, normed_kerr, 'black');
hold on
semilogy(e_pump_ranges, normed_brunel, 'b-.');
semilogy(e_pump_ranges, normed_injection, 'b-');