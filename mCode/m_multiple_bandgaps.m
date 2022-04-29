clear all;
c = 299792458;
q = 1.60217662e-19;
me = 9.10938e-31;
n0 = 2.2e28; 
bandgaps = [7.5, 7.7, 8];
wavelength_probe = 800e-9;
wavelength_pump = 2100e-9;

% integration params
t_end = 1000e-15;
delta_t = 5e-18;
t = 0:delta_t:t_end;
L = length(t);

%delay times
tau_pump = 430e-15;
tau_probe = 430e-15;


%varying pump intensities 
no_simulations = 35;
i_pump_ranges = linspace(3e16, 20e16, no_simulations);
amplitude_probe = intensity2amplitude(0.015e16); %0.015 TWcm^-2
amplitude_pump = intensity2amplitude(i_pump_ranges);
e_field_probe = zeros(3, L);
e_field_pump_parall = zeros(3, L);
e_field_pump_perpe = zeros(3, L);
normed_e_field_parall = zeros(1, L);
normed_e_field_perpe = zeros(1, L);


e_field_probe(1,:) = gaussian_efield_new(amplitude_probe, wavelength_probe, 45e-15, tau_probe, t);

injection_first_harm_parall = zeros(length(bandgaps),no_simulations);
brunel_first_harm_parall = zeros(length(bandgaps),no_simulations);
kerr_first_harm_parall = zeros(length(bandgaps),no_simulations);
overall_along_x_harm_parall = zeros(length(bandgaps),no_simulations);

injection_first_harm_perpe = zeros(length(bandgaps),no_simulations);
brunel_first_harm_perpe = zeros(length(bandgaps),no_simulations);
kerr_first_harm_perpe = zeros(length(bandgaps),no_simulations);
overall_along_x_harm_perpe = zeros(length(bandgaps),no_simulations);

third_term_parall = zeros(3, L);
third_term_perpe = zeros(3, L);

injection_first_harm_parall_ADK = zeros(length(bandgaps),no_simulations);
brunel_first_harm_parall_ADK = zeros(length(bandgaps),no_simulations);
kerr_first_harm_parall_ADK = zeros(length(bandgaps),no_simulations);
overall_along_x_harm_parall_ADK = zeros(length(bandgaps),no_simulations);

injection_first_harm_perpe_ADK = zeros(length(bandgaps),no_simulations);
brunel_first_harm_perpe_ADK = zeros(length(bandgaps),no_simulations);
kerr_first_harm_perpe_ADK = zeros(length(bandgaps),no_simulations);
overall_along_x_harm_perpe_ADK = zeros(length(bandgaps),no_simulations);

third_term_parall_ADK = zeros(3, L);
third_term_perpe_ADK = zeros(3, L);

n_fft = 2^nextpow2(L); %zero padding
f = 1/delta_t*(0:(n_fft/2))/n_fft;
delta_f = 1/delta_t;
f_pump = c / 2.1e-06;
f_probe = c / 8e-07;
first_harm = 2*f_pump + f_probe;
[~, idx] = min(abs(first_harm - f));
[~, idx_pump] = min(abs(f_pump - f));
[~, idx_probe] = min(abs(f_probe - f));


for b = 1:length(bandgaps)
    for i = 1:no_simulations
            e_field_pump_parall(1,:) = gaussian_efield_new(amplitude_pump(i), wavelength_pump, 140e-15, tau_pump, t);
            e_field_parall = e_field_pump_parall + e_field_probe;
            e_field_pump_perpe(2,:) = gaussian_efield_new(amplitude_pump(i), wavelength_pump, 140e-15, tau_pump, t);
            e_field_perpe = e_field_pump_perpe + e_field_probe;
            for j = 1:L
                normed_e_field_parall(:,j) = norm(e_field_parall(:,j));
                normed_e_field_perpe(:,j) = norm(e_field_perpe(:,j));
            end
            displacements_x_parall = displacement_x_new(bandgaps(b), normed_e_field_parall + 1, e_field_parall);
            ADK_parall = tangent_Gamma_ADK(normed_e_field_parall, bandgaps(b));
            rho_sfi_parall = integrate_population_cb(ADK_parall, delta_t, t);
            drho_parall = gradient(rho_sfi_parall, delta_t);
            third_term_parall(1,:) = gradient(displacements_x_parall(1,:).*drho_parall, delta_t);
            v0 = 0;
            brunel_current_density_parall = n0 * q * q/me * e_field_parall.*rho_sfi_parall;
            kerr_parall = kerr_current_binomial(e_field_pump_parall, e_field_probe, delta_t);
            injection_current_density_parall =  n0 .* q .* third_term_parall; 
            overall_current_x_parall = brunel_current_density_parall(1,:) + kerr_parall(1,:) + injection_current_density_parall(1,:);
    
            ft_brunel_current_parall = fft(brunel_current_density_parall(1,:), n_fft);
            ft_injection_current_parall = fft(injection_current_density_parall(1,:), n_fft);
            ft_kerr_current_parall = fft(kerr_parall(1,:), n_fft);
            ft_overall_current_x_parall = fft(overall_current_x_parall(1,:), n_fft);
            P_brunel_current_parall = abs(ft_brunel_current_parall/n_fft).^2;
            P_injection_current_parall = abs(ft_injection_current_parall/n_fft).^2;
            P_kerr_current_parall = abs(ft_kerr_current_parall/n_fft).^2;
            P_overall_current_parall = abs(ft_overall_current_x_parall/n_fft).^2;
            brunel_first_harm_parall(b, i) = P_brunel_current_parall(idx);
            kerr_first_harm_parall(b, i) = P_kerr_current_parall(idx);
            injection_first_harm_parall(b, i) = P_injection_current_parall(idx);
            overall_along_x_harm_parall(b, i) = P_overall_current_parall(idx);
            
            displacements_x_parall_ADK = displacement_x_new(bandgaps(b), normed_e_field_parall + 1, e_field_parall);
            ADK_parall_ADK = Gamma_ADK(bandgaps(b).*q, normed_e_field_parall, 1, 0, 0);
            rho_sfi_parall_ADK = integrate_population_cb(ADK_parall_ADK, delta_t, t);
            drho_parall_ADK = gradient(rho_sfi_parall_ADK, delta_t);
            third_term_parall_ADK(1,:) = gradient(displacements_x_parall_ADK(1,:).*drho_parall_ADK, delta_t);
            v0 = 0;
            brunel_current_density_parall_ADK = n0 * q * q/me * e_field_parall.*rho_sfi_parall_ADK;
            kerr_parall_ADK = kerr_current_binomial(e_field_pump_parall, e_field_probe, delta_t);
            injection_current_density_parall_ADK =  n0 .* q .* third_term_parall_ADK; 
            overall_current_x_parall_ADK = brunel_current_density_parall_ADK(1,:) + kerr_parall_ADK(1,:) + injection_current_density_parall_ADK(1,:);
    
            ft_brunel_current_parall_ADK = fft(brunel_current_density_parall_ADK(1,:), n_fft);
            ft_injection_current_parall_ADK = fft(injection_current_density_parall_ADK(1,:), n_fft);
            ft_kerr_current_parall_ADK = fft(kerr_parall_ADK(1,:), n_fft);
            ft_overall_current_x_parall_ADK = fft(overall_current_x_parall_ADK(1,:), n_fft);
            P_brunel_current_parall_ADK = abs(ft_brunel_current_parall_ADK/n_fft).^2;
            P_injection_current_parall_ADK = abs(ft_injection_current_parall_ADK/n_fft).^2;
            P_kerr_current_parall_ADK = abs(ft_kerr_current_parall_ADK/n_fft).^2;
            P_overall_current_parall_ADK = abs(ft_overall_current_x_parall_ADK/n_fft).^2;
            brunel_first_harm_parall_ADK(b, i) = P_brunel_current_parall_ADK(idx);
            kerr_first_harm_parall_ADK(b, i) = P_kerr_current_parall_ADK(idx);
            injection_first_harm_parall_ADK(b, i) = P_injection_current_parall_ADK(idx);
            overall_along_x_harm_parall_ADK(b, i) = P_overall_current_parall_ADK(idx);

            displacements_x_perpe = displacement_x_new(bandgaps(b), normed_e_field_perpe + 1, e_field_perpe);
            ADK_perpe = tangent_Gamma_ADK(normed_e_field_perpe, bandgaps(b));
            rho_sfi_perpe = integrate_population_cb(ADK_perpe, delta_t, t);
            drho_perpe = gradient(rho_sfi_perpe, delta_t);
            third_term_perpe(1,:) = gradient(displacements_x_perpe(1,:).*drho_perpe, delta_t);
            v0 = 0;
            brunel_current_density_perpe = n0 * q * q/me * e_field_perpe.*rho_sfi_perpe;
            kerr_perpe = kerr_current_binomial(e_field_pump_perpe, e_field_probe, delta_t);
            injection_current_density_perpe =  n0 .* q .* third_term_perpe; 
            overall_current_x_perpe = brunel_current_density_perpe(1,:) + kerr_perpe(1,:) + injection_current_density_perpe(1,:);
    
            ft_brunel_current_perpe = fft(brunel_current_density_perpe(1,:), n_fft);
            ft_injection_current_perpe = fft(injection_current_density_perpe(1,:), n_fft);
            ft_kerr_current_perpe = fft(kerr_perpe(1,:), n_fft);
            ft_overall_current_x_perpe = fft(overall_current_x_perpe(1,:), n_fft);
            P_brunel_current_perpe = abs(ft_brunel_current_perpe/n_fft).^2;
            P_injection_current_perpe = abs(ft_injection_current_perpe/n_fft).^2;
            P_kerr_current_perpe = abs(ft_kerr_current_perpe/n_fft).^2;
            P_overall_current_perpe = abs(ft_overall_current_x_perpe/n_fft).^2;
            brunel_first_harm_perpe(b, i) = P_brunel_current_perpe(idx);
            kerr_first_harm_perpe(b, i) = P_kerr_current_perpe(idx);
            injection_first_harm_perpe(b, i) = P_injection_current_perpe(idx);
            overall_along_x_harm_perpe(b, i) = P_overall_current_perpe(idx);
            
            displacements_x_perpe_ADK = displacement_x_new(bandgaps(b), normed_e_field_perpe + 1, e_field_perpe);
            ADK_perpe_ADK = Gamma_ADK(bandgaps(b).*q, normed_e_field_perpe, 1, 0, 0);
            rho_sfi_perpe_ADK = integrate_population_cb(ADK_perpe_ADK, delta_t, t);
            drho_perpe_ADK = gradient(rho_sfi_perpe_ADK, delta_t);
            third_term_perpe_ADK(1,:) = gradient(displacements_x_perpe_ADK(1,:).*drho_perpe_ADK, delta_t);
            v0 = 0;
            brunel_current_density_perpe_ADK = n0 * q * q/me * e_field_parall.*rho_sfi_perpe_ADK;
            kerr_perpe_ADK = kerr_current_binomial(e_field_pump_perpe, e_field_probe, delta_t);
            injection_current_density_perpe_ADK =  n0 .* q .* third_term_perpe_ADK; 
            overall_current_x_perpe_ADK = brunel_current_density_perpe_ADK(1,:) + kerr_perpe_ADK(1,:) + injection_current_density_perpe_ADK(1,:);
    
            ft_brunel_current_perpe_ADK = fft(brunel_current_density_perpe_ADK(1,:), n_fft);
            ft_injection_current_perpe_ADK = fft(injection_current_density_perpe_ADK(1,:), n_fft);
            ft_kerr_current_perpe_ADK = fft(kerr_perpe_ADK(1,:), n_fft);
            ft_overall_current_x_perpe_ADK = fft(overall_current_x_perpe_ADK(1,:), n_fft);
            P_brunel_current_perpe_ADK = abs(ft_brunel_current_perpe_ADK/n_fft).^2;
            P_injection_current_perpe_ADK = abs(ft_injection_current_perpe_ADK/n_fft).^2;
            P_kerr_current_perpe_ADK = abs(ft_kerr_current_perpe_ADK/n_fft).^2;
            P_overall_current_perpe_ADK = abs(ft_overall_current_x_perpe_ADK/n_fft).^2;
            brunel_first_harm_perpe_ADK(b, i) = P_brunel_current_perpe_ADK(idx);
            kerr_first_harm_perpe_ADK(b, i) = P_kerr_current_perpe_ADK(idx);
            injection_first_harm_perpe_ADK(b, i) = P_injection_current_perpe_ADK(idx);
            overall_along_x_harm_perpe_ADK(b, i) = P_overall_current_perpe_ADK(idx);
    end
end

figure(1);
% Get the handle of figure(n).
fig1_comps.fig = gcf;
hold on 
styles = {'-', '--', '-.'};
for b = 1:length(bandgaps)
    m = sqrt(overall_along_x_harm_parall(b, :)./overall_along_x_harm_perpe(b, :));
    p1 = plot(i_pump_ranges, m);
    p2 = plot(i_pump_ranges, sqrt(kerr_first_harm_parall(b, :)./kerr_first_harm_perpe(b, :)));
    p3 = plot(i_pump_ranges, sqrt((kerr_first_harm_parall(b, :) + brunel_first_harm_parall(b, :))./(kerr_first_harm_perpe(b, :) + brunel_first_harm_perpe(b, :))));
    set(p1, 'LineWidth', 1.5, 'color', 'blue', 'LineStyle', styles{b});
    set(p2, 'LineWidth', 1.5, 'color', 'black', 'LineStyle', styles{b});
    set(p3, 'LineWidth', 1.5, 'color', 'red', 'LineStyle', styles{b});
    m_ADK = sqrt(overall_along_x_harm_parall_ADK(b, :)./overall_along_x_harm_perpe_ADK(b, :));
    p4 = plot(i_pump_ranges, m_ADK);
    p5 = plot(i_pump_ranges, sqrt((kerr_first_harm_parall_ADK(b, :) + brunel_first_harm_parall_ADK(b, :))./(kerr_first_harm_perpe_ADK(b, :) + brunel_first_harm_perpe_ADK(b, :))));
    set(p4, 'LineWidth', 0.5, 'color', 'black', 'LineStyle', styles{b});
    set(p5, 'LineWidth', 0.5, 'color', 'black', 'LineStyle', styles{b});
end
fig1_comps.tile1.plotXLabel = xlabel('$$I_{pump}$$ in $$Wm^{-2}$$');
fig1_comps.tile1.plotYLabel = ylabel('$$m = \sqrt{I_{\mid\mid}/I_{\perp}}$$');
%title('$$E_{Gap}=7.5 eV , \Delta t = 1e-18s, \tau_{delay} = 0$$', 'Interpreter','latex')
legend('m', 'Kerr', 'Kerr+Brunel', 'location', 'northwest')
ylim([2,14])
grid on
xlim([i_pump_ranges(1), i_pump_ranges(end)])
STANDARDIZE_FIGURE(fig1_comps)