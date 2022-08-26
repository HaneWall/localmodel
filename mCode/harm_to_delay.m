%script that does compare the first harmonic n = 1 for 2*n omega_pump +
%omega_probe in respect to the delay (parallel polarization)
clear all;
%physical constants
c = 299792458;
q = 1.60217662e-19;
me = 9.10938e-31;
n0 = 2.2e28;                                     %molecular density for si02

phase_plot=false;
hue_plot=false;
delay_plot_all=true;
delay_filter_plot=false;
individual_phase_plot = false;
sawtooth_plot_phases = false;
finer_sawtooth_plot_phases = false;
stft_plot = false;
cwt_plot = false;
hilbert_phase_plot = false;

%simulation parameters 
bandgap = 7.5; % in eV

amplitude_probe= intensity2amplitude(0.015e16);  %0.015 TWcm^-2
amplitude_pump = intensity2amplitude(12e16);     %12 TWcm^-2

wavelength_probe = 800e-9;
wavelength_pump = 2100e-9;

% FWHM in terms of intensity 
fwhm_pump = 140e-15;
fwhm_probe = 45e-15;

%integration parameters
t_end = 1000e-15;
delta_t = 3.8e-18;
t = 0:delta_t:t_end;

%delay times
L_Delay = 51;
delay_between_pulses = linspace(-200e-15, 200e-15, L_Delay); %80
pump_peak_t = 500e-15;
tau_pump = ones(1, L_Delay) * pump_peak_t;
tau_probe = tau_pump + delay_between_pulses;

%allocate memory for first harmonic
L = length(t);
n_fft = 2^nextpow2(L); %zero padding
power_density_harmonic = zeros(length(delay_between_pulses), n_fft/2 + 1);
power_kerr_harmonic = zeros(length(delay_between_pulses), n_fft/2 + 1);
power_injection_harmonic = zeros(length(delay_between_pulses), n_fft/2 + 1);
power_brunel_harmonic = zeros(length(delay_between_pulses), n_fft/2 + 1);

phase_density_harmonic = zeros(length(delay_between_pulses), n_fft/2 + 1);
phase_kerr_harmonic = zeros(length(delay_between_pulses), n_fft/2 + 1);
phase_injection_harmonic = zeros(length(delay_between_pulses), n_fft/2 + 1);
phase_brunel_harmonic = zeros(length(delay_between_pulses), n_fft/2 + 1);

kerr_container = zeros(length(delay_between_pulses), n_fft);
injection_container = zeros(length(delay_between_pulses), n_fft);
brunel_container = zeros(length(delay_between_pulses), n_fft);
overall_container = zeros(length(delay_between_pulses), n_fft);

ft_overall_current_pump_probe = zeros(length(delay_between_pulses), n_fft/2 + 1);
overall_current_pump_probe = zeros(length(delay_between_pulses), L);

e_field_probe = zeros(3, L);
e_field_pump = zeros(3, L);
third_term = zeros(3, L);
normed_e_field = zeros(1, L);
t_all = 0:delta_t:(n_fft-1)*delta_t;

for i = 1:length(delay_between_pulses)
    e_field_pump(1, :) = gaussian_efield_new(amplitude_pump, wavelength_pump, fwhm_pump, tau_pump(i), t);
    e_field_probe(1, :) = gaussian_efield_new(amplitude_probe, wavelength_probe, fwhm_probe, tau_probe(i), t);
    e_field = e_field_pump + e_field_probe;
    for j = 1:L
        normed_e_field(:,j) = norm(e_field(:,j));
    end
    displacements_x = displacement_x_new(bandgap, abs(normed_e_field) + 1, e_field);
    ADK = tangent_Gamma_ADK(normed_e_field, bandgap);
    
    rho_sfi = integrate_population_cb(ADK, delta_t, t);
    drho = gradient(rho_sfi, delta_t);
    third_term(1, :) = gradient(displacements_x(1, :).*drho, delta_t);
    v0 = 0;
    plasma_current_density = n0 .* q .* (q/me .* e_field.*rho_sfi + v0*drho + third_term);
    brunel_current_density = n0 * q * q/me * e_field.*rho_sfi;
    kerr = kerr_current_binomial(e_field_probe, e_field_pump, delta_t);  
    injection_current_density =  n0 .* q .* third_term; 
    overall_current_density = brunel_current_density(1,:) + kerr(1,:) + injection_current_density(1,:);
    
    ft_plasma_current = fft(plasma_current_density(1, :),n_fft);
    ft_brunel_current = fft(brunel_current_density(1, :), n_fft);
    ft_injection_current = fft(injection_current_density(1, :), n_fft);
    ft_kerr_current = fft(kerr(1, :), n_fft); 
    ft_overall_current = fft(overall_current_density(1, :), n_fft);
    ft_overall_current_pump_probe(i, :) = ft_overall_current(1:n_fft/2 + 1);
    overall_current_pump_probe(i, :) = overall_current_density(1, :);
    
    kerr_power_spec = abs(ft_kerr_current/n_fft).^2;
    injection_power_spec = abs(ft_injection_current/n_fft).^2;
    whole_power_spec = abs(ft_overall_current/n_fft).^2;
    brunel_power_spec = abs(ft_brunel_current/n_fft).^2;
    
    brunel_phase_spec = angle(ft_brunel_current/n_fft);
    injection_phase_spec = angle(ft_injection_current/n_fft);
    kerr_phase_spec = angle(ft_kerr_current/n_fft);
    overall_phase_spec = angle(ft_overall_current/n_fft);

    power_injection_harmonic(i,:) = injection_power_spec(1:n_fft/2 + 1);
    power_kerr_harmonic(i,:) = kerr_power_spec(1:n_fft/2 + 1);
    power_density_harmonic(i,:) = whole_power_spec(1:n_fft/2 + 1);
    power_brunel_harmonic(i,:) = brunel_power_spec(1:n_fft/2 + 1);

    phase_injection_harmonic(i,:) = injection_phase_spec(1:n_fft/2 + 1);
    phase_kerr_harmonic(i,:) = kerr_phase_spec(1:n_fft/2 + 1);
    phase_density_harmonic(i,:) = overall_phase_spec(1:n_fft/2 + 1);
    phase_brunel_harmonic(i,:) = brunel_phase_spec(1:n_fft/2 + 1);

    kerr_container(i, 1:L) = kerr(1,:);
    injection_container(i, 1:L) = injection_current_density(1,:);
    brunel_container(i, 1:L) = brunel_current_density(1, :);
    overall_container(i, 1:L) = overall_current_density(1, :);
end

log_power_spec = log10(power_density_harmonic);
log_kerr_spec = log10(power_kerr_harmonic);
log_injection_spec = log10(power_injection_harmonic);
log_brunel_spec = log10(power_brunel_harmonic);
f = 1/delta_t*(0:(n_fft/2))/n_fft;
delta_f = 1/delta_t;
f_pump = c / 2.1e-06;
f_probe = c / 8e-07;
first_harm = 2*f_pump + f_probe;
[~, idx] = min(abs(first_harm - f));
[~, idx_pump] = min(abs(f_pump - f));
[~, idx_probe] = min(abs(f_probe - f));

if delay_filter_plot
   %modify log power spec at each harmonic of omega pump 
    %idx_half_band = 55;
    
    power_density_harmonic_pump = zeros(length(delay_between_pulses), n_fft/2 + 1);
    power_kerr_harmonic_pump = zeros(length(delay_between_pulses), n_fft/2 + 1);
    power_injection_harmonic_pump = zeros(length(delay_between_pulses), n_fft/2 + 1);
    power_brunel_harmonic_pump = zeros(length(delay_between_pulses), n_fft/2 + 1);

    ft_overall_current_pump = zeros(length(delay_between_pulses), n_fft/2 + 1);
    ft_overall_current_diff_pump = zeros(length(delay_between_pulses), n_fft/2 + 1);
    diff_overall_current = zeros(length(delay_between_pulses), L);

    for i = 1:length(delay_between_pulses)
        e_field_pump(1, :) = gaussian_efield_new(amplitude_pump, wavelength_pump, fwhm_pump, tau_pump(i), t);
        e_field_probe(1, :) = zeros(1, length(t));
        e_field = e_field_pump + e_field_probe;
        for j = 1:L
            normed_e_field(:,j) = norm(e_field(:,j));
        end
        displacements_x = displacement_x_new(bandgap, abs(normed_e_field) + 1, e_field);
        ADK = tangent_Gamma_ADK(normed_e_field, bandgap);
    
        rho_sfi = integrate_population_cb(ADK, delta_t, t);
        drho = gradient(rho_sfi, delta_t);
        third_term(1, :) = gradient(displacements_x(1, :).*drho, delta_t);
        v0 = 0;
        plasma_current_density = n0 * q * (q/me * e_field.*rho_sfi + v0*drho + third_term);
        brunel_current_density = n0 * q * q/me * e_field.*rho_sfi;
        kerr = kerr_current_binomial(e_field_pump, e_field_probe, delta_t);  
        injection_current_density =  n0 * q * third_term; 
        overall_current_density = brunel_current_density(1, :) + kerr(1, :) + injection_current_density(1, :);
        
        ft_plasma_current = fft(plasma_current_density(1, :),n_fft);
        ft_brunel_current = fft(brunel_current_density(1, :), n_fft);
        ft_injection_current = fft(injection_current_density(1, :), n_fft);
        ft_kerr_current = fft(kerr(1, :), n_fft); 
        
        ft_overall_current = fft(overall_current_density, n_fft);
        ft_overall_current_pump(i, :) = ft_overall_current(1:n_fft/2 + 1);

        diff_overall_current(i, :) = overall_current_pump_probe(i, :) - overall_current_density(1,:);
        
        ft_overall_current_diff = fft(diff_overall_current(i, :), n_fft);
        ft_overall_current_diff_pump(i,:) = ft_overall_current_diff(1:n_fft/2 + 1);
        
        kerr_power_spec = abs(ft_kerr_current/n_fft).^2;
        injection_power_spec = abs(ft_injection_current/n_fft).^2;
        whole_power_spec = abs(ft_overall_current/n_fft).^2;
        brunel_power_spec = abs(ft_brunel_current/n_fft).^2;
    
        power_injection_harmonic_pump(i,:) = injection_power_spec(1:n_fft/2 + 1);
        power_kerr_harmonic_pump(i,:) = kerr_power_spec(1:n_fft/2 + 1);
        power_density_harmonic_pump(i,:) = whole_power_spec(1:n_fft/2 + 1);
        power_brunel_harmonic_pump(i,:) = brunel_power_spec(1:n_fft/2 + 1);

    end
    
    %for i=0:7
       % modified_log_spec(:, (2*i+1)*idx_pump - idx_half_band:(2*i+1)*idx_pump + idx_half_band) = 30; 
       %modified_log_spec(:, (2*i+1)*idx_probe - idx_half_band:(2*i+1)*idx_probe + idx_half_band) = 30;
    %end
    %min = abs(min(min(power_density_harmonic - power_brunel_harmonic_pump)));
    
    figure(6)
    % Get the handle of figure(n).
    %modified_log_spec = log10(abs(ft_overall_current_diff_pump/n_fft).^2);
    modified_log_spec_img = log10(abs(abs(ft_overall_current_pump_probe/n_fft).^2 - abs(ft_overall_current_pump/n_fft).^2));
    fig1_comps.fig = gcf;
    p1 = imagesc(f./f_pump, delay_between_pulses.*10^(15), modified_log_spec_img, [47.5, 52]);
    fig1_comps.tile1.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig1_comps.tile1.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), modified_log_spec_img, [49 50 51]);
    c1.LineColor = 'black';
    c1.LineWidth = 0.5;
    caxis(cRange)
    title('Overall')
    xline(f(idx_pump)/f_pump, 'w-.');
    xline(f(idx_probe)/f_pump, 'w--');
    for i=1:5
        xline(f(2*i*idx_pump + idx_probe)/f_pump, 'w-');
    end
    xlim([0, 11] );
    ax = gca;
    ax.YDir = 'normal';
    colormap jet;
    STANDARDIZE_FIGURE(fig1_comps);
end


if delay_plot_all
%     figure(7)
%     % Get the handle of figure(n).
%     fig2_comps.fig = gcf;
%     imagesc(f./f_pump, delay_between_pulses.*10^(15), log_power_spec, [47.5, 52]);
%     fig2_comps.tile1.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
%     fig2_comps.tile1.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
%     title('$$\log_{10}|\mathcal{F}(\partial_t \mathbf{j}})|^2$$')
%     xline(f(idx)/f_pump, 'w-');
%     xline(f(idx_pump)/f_pump, 'w-.');
%     xline(f(idx_probe)/f_pump, 'w--');
%     xlim([0, 11]);
%     ax = gca;
%     ax.YDir = 'normal';
%     colormap jet;
%     STANDARDIZE_FIGURE(fig2_comps);
    
    figure(8);
    fig3_comps.fig = gcf;
    fig3_comps.t1 = tiledlayout(fig3_comps.fig, 1, 4);
    fig3_comps.n(1) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), log_power_spec, [47.5, 52]);
    colormap jet;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'black';
    c1.LineWidth = 0.5;
    caxis(cRange)
    title('$$\partial_t j_{Overall}(\omega)$$')
    xline(f(idx)/f_pump, 'w-');
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    ax = gca;
    ax.YDir = 'normal';
    fig3_comps.tile1.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig3_comps.tile1.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');

    

    %STANDARDIZE_FIGURE(fig3_comps);

    fig3_comps.n(2) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), log_injection_spec, [47.5, 52]);
    colormap jet;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_injection_spec, [49 50 51]);
    c1.LineColor = 'black';
    c1.LineWidth = 0.5;
    caxis(cRange)
    title('$$\partial_t j_{Injection}(\omega)$$')
    xline(f(idx)/f_pump, 'w-');
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    ax = gca;
    ax.YDir = 'normal';
    fig3_comps.tile2.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    
    fig3_comps.n(3) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), log_kerr_spec, [47.5, 52]);
    colormap jet;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_kerr_spec, [49 50 51]);
    c1.LineColor = 'black';
    c1.LineWidth = 0.5;
    caxis(cRange)
    title('$$\partial_t j_{Kerr}(\omega)$$')
    xline(f(idx)/f_pump, 'w-');
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    ax = gca;
    ax.YDir = 'normal';
    fig3_comps.tile3.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    
    
    fig3_comps.n(4) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), log_brunel_spec, [47.5, 52]);
    colormap jet;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_brunel_spec, [49 50 51]);
    c1.LineColor = 'black';
    c1.LineWidth = 0.5;
    caxis(cRange)
    title('$$\partial_t j_{Brunel}(\omega)$$')
    xline(f(idx)/f_pump, 'w-');
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    ax = gca;
    ax.YDir = 'normal';
    fig3_comps.tile4.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    cbh = colorbar;
    cbh.Layout.Tile = 'north'; 

    STANDARDIZE_FIGURE(fig3_comps);
end

if phase_plot
    figure(3);

    fig4_comps.fig = gcf;
    fig4_comps.t1 = tiledlayout(fig4_comps.fig, 1, 3);

    fig4_comps.n(1) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), mod(phase_injection_harmonic - phase_brunel_harmonic + pi,2*pi) - pi)
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    title('$$\angle\partial_t j_{Injection} - \angle\partial_t j_{Brunel}$$');
    ax = gca;
    ax.YDir = 'normal';
    fig4_comps.tile1.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig4_comps.tile1.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
    
    fig4_comps.n(2) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), mod(phase_injection_harmonic - phase_kerr_harmonic + pi,2*pi) - pi)
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    [M2,c2] = contour(f/f_pump, delay_between_pulses.*10^(15), log_kerr_spec, [49 50 51]);
    c2.LineColor = 'black';
    c2.LineWidth = 1;
    [M3,c3] = contour(f/f_pump, delay_between_pulses.*10^(15), log_injection_spec, [49 50 51]);
    c3.LineColor = 'blue';
    c3.LineWidth = 1;
    caxis(cRange)
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    title('$$\angle\partial_t j_{Injection} - \angle\partial_t j_{Kerr}$$');
    ax = gca;
    ax.YDir = 'normal';
    fig4_comps.tile2.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig4_comps.tile2.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');

    fig4_comps.n(3) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), mod(phase_brunel_harmonic - phase_kerr_harmonic + pi,2*pi) - pi)
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    title('$$\angle\partial_t j_{Brunel} - \angle\partial_t j_{Kerr}$$');
    ax = gca;
    ax.YDir = 'normal';
    fig4_comps.tile3.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig4_comps.tile3.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
    cbh = colorbar;
    cbh.Ticks = linspace(-pi, pi, 9);
    cbh.TickLabels = {'$-\pi$', '$- 3/4 \pi$', '$-\pi/2$', '$-\pi/4$', '$0$', '$\pi/4$', '$\pi/2$', '$3/4\pi$', '$\pi$'};
    cbh.TickLabelInterpreter = 'latex';
    STANDARDIZE_FIGURE(fig4_comps);

    figure(5);
    fig5_comps.fig = gcf;
    fig5_comps.t1 = tiledlayout(fig5_comps.fig, 1, 3);
    
    fig5_comps.n(1) = nexttile;

    imagesc(f/f_pump, delay_between_pulses.*10^(15), mod(phase_density_harmonic - phase_kerr_harmonic + pi,2*pi) - pi)
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    title('Phase zw. Overall & Kerr')
    ax = gca;
    ax.YDir = 'normal';
    fig5_comps.tile1.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig5_comps.tile1.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
    
    fig5_comps.n(2) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), mod(phase_density_harmonic - phase_brunel_harmonic + pi,2*pi) - pi)
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    title('Phase zw. Overall & Brunel')
    ax = gca;
    ax.YDir = 'normal';
    fig5_comps.tile2.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig5_comps.tile2.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');

    fig5_comps.n(3) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), mod(phase_density_harmonic - phase_injection_harmonic + pi,2*pi) - pi)
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    title('Phase zw. Overall & Injection')
    ax = gca;
    ax.YDir = 'normal';
    fig5_comps.tile3.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig5_comps.tile3.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
    cbh = colorbar;
    cbh.Ticks = linspace(-pi, pi, 9);
    cbh.TickLabels = {'$-\pi$', '$- 3/4 \pi$', '$-\pi/2$', '$-\pi/4$', '$0$', '$\pi/4$', '$\pi/2$', '$3/4\pi$', '$\pi$'};
    cbh.TickLabelInterpreter = 'latex';
    STANDARDIZE_FIGURE(fig5_comps);

end

if sawtooth_plot_phases
    figure(8)
    fig8_comps.fig = gcf;
    fig8_comps.t1  = tiledlayout(fig8_comps.fig, 1, 2);
    fig8_comps.n(1) = nexttile;
    title('$$\angle\partial_t j_{Injection} - \angle\partial_t j_{Kerr}$$')
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    hold on 
    for i = 11:20:101
       plot(f/f_pump, mod(phase_injection_harmonic(i, :) - phase_kerr_harmonic(i, :) + pi, 2*pi) - pi, "LineWidth", 1)
    end
    fig8_comps.tile1.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    
    xline(f(idx)/f_pump, 'b-');

    fig8_comps.n(2) = nexttile;
    title('$$\angle\partial_t j_{Brunel} - \angle\partial_t j_{Kerr}$$')
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    hold on 
    for i = 11:20:101
       plot(f/f_pump, mod(phase_brunel_harmonic(i, :) - phase_kerr_harmonic(i, :) + pi, 2*pi) - pi, "LineWidth", 1)
    end
    fig8_comps.tile2.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    xline(f(idx)/f_pump, 'b-');
    leg = legend("$$\tau_{Delay}=$$" + delay_between_pulses(11)*10^(15)+"fs", "$$\tau_{Delay}=$$" + delay_between_pulses(31)*10^(15)+"fs", "$$\tau_{Delay}=$$" + delay_between_pulses(51)*10^(15) +"fs", "$$\tau_{Delay}=$$" + delay_between_pulses(71)*10^(15)+"fs", "$$\tau_{Delay}=$$" + delay_between_pulses(91)*10^(15)+"fs");
    leg.Layout.Tile = 'north'; 
    STANDARDIZE_FIGURE(fig8_comps)
end

if finer_sawtooth_plot_phases
    figure(8)
    fig8_comps.fig = gcf;
    fig8_comps.t1  = tiledlayout(fig8_comps.fig, 1, 2);
    fig8_comps.n(1) = nexttile;
    title('$$\angle\partial_t j_{Injection} - \angle\partial_t j_{Kerr}$$')
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    hold on 
    for i = 51:10:91
       plot(f/f_pump, mod(phase_injection_harmonic(i, :) - phase_kerr_harmonic(i, :) + pi, 2*pi) - pi, "LineWidth", 1)
    end
    fig8_comps.tile1.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    
    xline(f(idx)/f_pump, 'b-');

    fig8_comps.n(2) = nexttile;
    title('$$\angle\partial_t j_{Brunel} - \angle\partial_t j_{Kerr}$$')
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    hold on 
    for i = 51:10:91
       plot(f/f_pump, mod(phase_brunel_harmonic(i, :) - phase_kerr_harmonic(i, :) + pi, 2*pi) - pi, "LineWidth", 1)
    end
    fig8_comps.tile2.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    xline(f(idx)/f_pump, 'b-');
    leg = legend("$$\tau_{Delay}=$$" + delay_between_pulses(51)*10^(15)+"fs", "$$\tau_{Delay}=$$" + delay_between_pulses(61)*10^(15)+"fs", "$$\tau_{Delay}=$$" + delay_between_pulses(71)*10^(15) +"fs", "$$\tau_{Delay}=$$" + delay_between_pulses(81)*10^(15)+"fs", "$$\tau_{Delay}=$$" + delay_between_pulses(91)*10^(15)+"fs");
    leg.Layout.Tile = 'north'; 
    STANDARDIZE_FIGURE(fig8_comps)
end


if hilbert_phase_plot
    phase_hilbert_kerr = zeros(length(delay_between_pulses), n_fft);
    phase_hilbert_brunel = zeros(length(delay_between_pulses), n_fft);
    phase_hilbert_injection = zeros(length(delay_between_pulses), n_fft);
    phase_hilbert_overall = zeros(length(delay_between_pulses), n_fft);

    for N = 1:1:length(delay_between_pulses)
        [~, phase_hilbert_kerr(N, :)] = analytic_signal_with_t(kerr_container(N, :), t_all, pump_peak_t);
        [~, phase_hilbert_brunel(N, :)] = analytic_signal_with_t(brunel_container(N, :), t_all, pump_peak_t);
        [~, phase_hilbert_injection(N, :)] = analytic_signal_with_t(injection_container(N, :), t_all, pump_peak_t);
        [~, phase_hilbert_overall(N, :)] = analytic_signal_with_t(overall_container(N, :), t_all, pump_peak_t);
    end

    
%     [~, phase_hilbert_kerr(1:N, :)] = analytic_signal_with_t(kerr_container(1:N, :), t_all, pump_peak_t);
%     [~, phase_hilbert_brunel(1:N, :)] = analytic_signal_with_t(brunel_container(1:N, :), t_all, pump_peak_t);
%     [~, phase_hilbert_injection(1:N, :)] = analytic_signal_with_t(injection_container(1:N, :), t_all, pump_peak_t);
%     [~, phase_hilbert_overall(1:N, :)] = analytic_signal_with_t(overall_container(1:N, :), t_all, pump_peak_t);
    
    f_shifted = 1/(n_fft*delta_t)*(-(floor((n_fft)/2)):1:(floor((n_fft-1)/2)));
    middle = ceil(n_fft/2);
    [~, first_idx_s] = min(abs(first_harm + f_shifted(1:middle)));
    second_idx_s = first_idx_s + 2*(middle-first_idx_s);

    figure(7)
    fig7_comps.fig = gcf;
    fig7_comps.t1 = tiledlayout(fig7_comps.fig, 1, 4);

    fig7_comps.n(1) = nexttile;
    imagesc(f_shifted/f_pump, delay_between_pulses.*10^(15), phase_hilbert_overall);
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f_shifted(middle:end)/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f_shifted(second_idx_s - 60)/f_pump, f_shifted(second_idx_s + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    ax = gca;
    ax.YDir = 'normal';
    title('$$\angle\partial_t j$$')
    fig7_comps.tile1.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
    fig7_comps.tile1.plotXLabel = xlabel('$$\omega/\omega_{pump}$$');
    

    fig7_comps.n(2) = nexttile;
    imagesc(f_shifted/f_pump, delay_between_pulses.*10^(15), phase_hilbert_injection);
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f_shifted(middle:end)/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f_shifted(second_idx_s - 60)/f_pump, f_shifted(second_idx_s + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    ax = gca;
    ax.YDir = 'normal';
    title('$$\angle\partial_t j_{Injection}$$')
 

    fig7_comps.n(3) = nexttile;
    imagesc(f_shifted/f_pump, delay_between_pulses.*10^(15), phase_hilbert_kerr);
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f_shifted(middle:end)/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f_shifted(second_idx_s - 60)/f_pump, f_shifted(second_idx_s + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    ax = gca;
    ax.YDir = 'normal';
    title('$$\angle\partial_t j_{Kerr}$$')
    fig7_comps.tile3.plotXLabel = xlabel('$$\omega/\omega_{pump}$$');

    fig7_comps.n(4) = nexttile;
    imagesc(f_shifted/f_pump, delay_between_pulses.*10^(15), phase_hilbert_brunel);
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f_shifted(middle:end)/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f_shifted(second_idx_s - 60)/f_pump, f_shifted(second_idx_s + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    ax = gca;
    ax.YDir = 'normal';
    title('$$\angle\partial_t j_{Brunel}$$')
    fig7_comps.tile4.plotXLabel = xlabel('$$\omega/\omega_{pump}$$');

    cbh = colorbar;
    cbh.Ticks = linspace(-pi, pi, 9);
    cbh.TickLabels = {'$-\pi$', '$- 3/4 \pi$', '$-\pi/2$', '$-\pi/4$', '$0$', '$\pi/4$', '$\pi/2$', '$3/4\pi$', '$\pi$'};
    cbh.TickLabelInterpreter = 'latex';
    cbh.Layout.Tile = 'south'; 

    STANDARDIZE_FIGURE(fig7_comps);
end

if individual_phase_plot
    figure(6);
    fig6_comps.fig = gcf;
    fig6_comps.t1 = tiledlayout(fig6_comps.fig, 2, 2);
    fig6_comps.n(1) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), mod(phase_density_harmonic + pi,2*pi) - pi);
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    title('$$\angle\partial_t j$$')
    ax = gca;
    ax.YDir = 'normal';
    fig6_comps.tile1.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');

    fig6_comps.n(2) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), mod(phase_kerr_harmonic + pi,2*pi) - pi);
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    title('$$\angle\partial_t j_{Kerr}$$')
    ax = gca;
    ax.YDir = 'normal';
    %fig6_comps.tile2.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');

    fig6_comps.n(3) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), mod(phase_injection_harmonic + pi,2*pi) - pi);
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    title('$$\angle\partial_t j_{Injection}$$')
    ax = gca;
    ax.YDir = 'normal';
    fig6_comps.tile3.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
    fig6_comps.tile3.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');

    fig6_comps.n(4) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), mod(phase_brunel_harmonic + pi,2*pi) - pi);
    colormap hsv;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_power_spec, [49 50 51]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    title('$$\angle\partial_t j_{Brunel}$$')
    ax = gca;
    ax.YDir = 'normal';
    fig6_comps.tile4.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');

    cbh = colorbar;
    cbh.Ticks = linspace(-pi, pi, 9);
    cbh.TickLabels = {'$-\pi$', '$- 3/4 \pi$', '$-\pi/2$', '$-\pi/4$', '$0$', '$\pi/4$', '$\pi/2$', '$3/4\pi$', '$\pi$'};
    cbh.TickLabelInterpreter = 'latex';
    cbh.Layout.Tile = 'east'; 
    STANDARDIZE_FIGURE(fig6_comps)
end

if stft_plot
    N = length(delay_between_pulses);
    delay_idx = ceil(3*N/8);
    % define analysis parameters
    x_overall = overall_container(delay_idx, :);
    x_kerr = kerr_container(delay_idx, :);
    x_brunel = brunel_container(delay_idx, :);
    x_injection = injection_container(delay_idx, :);
    fs = delta_f;
    
    wlen = 2^13;                        % window length (recomended to be power of 2)
    hop = wlen/4;                       % hop size (recomended to be power of 2)
    %nfft = 4096;                        % number of fft points (recomended to be power of 2)
    
    % perform STFT
    win = balckman(wlen, 'periodic');
    [S_overall, f, t] = stft(x_overall, win, hop, n_fft, fs);
    [S_kerr, f, t] = stft(x_kerr, win, hop, n_fft, fs);
    [S_brunel, f, t] = stft(x_brunel, win, hop, n_fft, fs);
    [S_injection, f, t] = stft(x_injection, win, hop, n_fft, fs);
    
    % calculate the coherent amplification of the window
    C = sum(win)/wlen;
    
    % take the amplitude of fft(x) and scale it, so not to be a
    % function of the length of the window and its coherent amplification
    S_overall = abs(S_overall)/wlen/C;
    S_kerr = abs(S_kerr)/wlen/C;
    S_brunel = abs(S_brunel)/wlen/C;
    S_injection = abs(S_injection)/wlen/C;
    
    % correction of the DC & Nyquist component
    if rem(n_fft, 2)                     % odd nfft excludes Nyquist point
        S_overall(2:end, :) = S_overall(2:end, :).*2;
        S_kerr(2:end, :) = S_kerr(2:end, :).*2;
        S_brunel(2:end, :) = S_brunel(2:end, :).*2;
        S_injection(2:end, :) = S_injection(2:end, :).*2;
    else                                % even nfft includes Nyquist point
        S_overall(2:end-1, :) = S_overall(2:end-1, :).*2;
        S_kerr(2:end-1, :) = S_kerr(2:end-1, :).*2;
        S_injection(2:end-1, :) = S_injection(2:end-1, :).*2;
        S_brunel(2:end-1, :) = S_brunel(2:end-1, :).*2;
    end
    
    % convert amplitude spectrum to dB (min = -120 dB)
    S_overall = log10(S_overall + 1e-6);
    S_injection = log10(S_injection + 1e-6);
    S_brunel = log10(S_brunel + 1e-6);
    S_kerr = log10(S_kerr + 1e-6);

    % plot the spectrogram
    figure(10)
    fig10_comps.fig = gcf;
    fig10_comps.t1 = tiledlayout(fig10_comps.fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    fig10_comps.n(1) = nexttile;
    imagesc(t.*10^15, f./f_pump, S_overall)
    ax = gca;
    set(ax, 'YDir', 'normal')
    xlim([50, 950]) 
    yline(first_harm/f_pump, 'color', 'w', 'LineWidth', 1.5);
    ylim([0, 10])
    caxis([24.7, 27.7])
    colormap jet;
    %xlabel('$$t$$ in fs')
    ylabel('$$\omega / \omega_{pump}$$')
    title("$$\partial_t j_{overall}, \tau_{delay} = $$" + delay_between_pulses(delay_idx)*10^15 + "fs")
    
    fig10_comps.n(2) = nexttile;
    imagesc(t.*10^15, f./f_pump, S_kerr)
    ax = gca;
    set(ax, 'YDir', 'normal')
    xlim([50, 950]) 
    yline(first_harm/f_pump, 'color', 'w', 'LineWidth', 1.5);
    ylim([0, 10])
    caxis([24.7, 27.7])
    colormap jet;
    %xlabel('$$t$$ in fs')
    %ylabel('$$\omega / \omega_{pump}$$')
    title("$$\partial_t j_{kerr}, \tau_{delay} = $$" + delay_between_pulses(delay_idx)*10^15 + "fs")

    fig10_comps.n(3) = nexttile;
    imagesc(t.*10^15, f./f_pump, S_brunel)
    ax = gca;
    set(ax, 'YDir', 'normal')
    xlim([50, 950]) 
    yline(first_harm/f_pump, 'color', 'w', 'LineWidth', 1.5);
    ylim([0, 10])
    caxis([24.7, 27.7])
    colormap jet;
    xlabel('$$t$$ in fs')
    ylabel('$$\omega / \omega_{pump}$$')
    title("$$\partial_t j_{brunel}, \tau_{delay} = $$" + delay_between_pulses(delay_idx)*10^15 + "fs")

    fig10_comps.n(4) = nexttile;
    imagesc(t.*10^15, f./f_pump, S_injection)
    ax = gca;
    set(ax, 'YDir', 'normal')
    xlim([50, 950]) 
    yline(first_harm/f_pump, 'color', 'w', 'LineWidth', 1.5);
    ylim([0, 10])
    caxis([24.7, 27.7])
    colormap jet;
    xlabel('$$t$$ in fs')
    %ylabel('$$\omega / \omega_{pump}$$')
    title("$$\partial_t j_{injection}, \tau_{delay} = $$" + delay_between_pulses(delay_idx)*10^15 + "fs")
    STANDARDIZE_FIGURE(fig10_comps);
end


if cwt_plot
    delay_idx = 81;
    [ovr_cwt, freq] = cwt(overall_container(delay_idx, :), delta_f);

    kerr_cwt = cwt(kerr_container(delay_idx, :), delta_f);

    brunel_cwt = cwt(brunel_container(delay_idx, :), delta_f);

    injection_cwt = cwt(injection_container(delay_idx, :), delta_f);

    p_ovr_cwt = log10(abs(ovr_cwt));
    p_kerr_cwt = log10(abs(kerr_cwt));
    p_brunel_cwt = log10(abs(brunel_cwt));
    p_injection_cwt = log10(abs(injection_cwt));

    order = freq./f_pump;
    harmonic = (2*f_pump + f_probe)/f_pump;

    figure(10)
    fig10_comps.fig = gcf;
    fig10_comps.t1 = tiledlayout(fig10_comps.fig, 2, 2);
    fig10_comps.n(1) = nexttile;
    imagesc((0:delta_t:delta_t*n_fft)*10^15, order, p_ovr_cwt);
    colormap turbo;
    hold on 
    yline(harmonic, 'w-')
    xline(tau_pump(delay_idx)*10^15, 'w-')
    xline(tau_probe(delay_idx)*10^15, 'w:')
    caxis([23, 30]);
    title(['$$ \partial_t j, \tau_{Delay}=$$', num2str(delay_between_pulses(delay_idx)), ' fs'])
    xlim([0, 1000])
    ylabel('$$f/f_{pump}$$')
    ax = gca;
    set(ax, 'YDir', 'normal');
    set(ax, 'YScale', 'log')

    fig10_comps.n(2) = nexttile;
    imagesc((0:delta_t:delta_t*n_fft)*10^15, order, p_kerr_cwt);
    colormap turbo;
    hold on 
    yline(harmonic, 'w-')
    xline(tau_pump(delay_idx)*10^15, 'w-')
    xline(tau_probe(delay_idx)*10^15, 'w:')
    caxis([23, 30]);
    title(['$$ \partial_t j_{Kerr}, \tau_{Delay}=$$', num2str(delay_between_pulses(delay_idx)), ' fs'])
    xlim([0, 1000])
    ax = gca;
    set(ax, 'YDir', 'normal');
    set(ax, 'YScale', 'log')

    fig10_comps.n(3) = nexttile;
    imagesc((0:delta_t:delta_t*n_fft)*10^15, order, p_brunel_cwt);
    colormap turbo;
    hold on 
    yline(harmonic, 'w-')
    xline(tau_pump(delay_idx)*10^15, 'w-')
    xline(tau_probe(delay_idx)*10^15, 'w:')
    caxis([23, 30]);
    title(['$$ \partial_t j_{Brunel}, \tau_{Delay}=$$', num2str(delay_between_pulses(delay_idx)), ' fs'])
    xlim([0, 1000])
    ylabel('$$f/f_{pump}$$')
    ax = gca;
    xlabel('t in fs')
    set(ax, 'YDir', 'normal');
    set(ax, 'YScale', 'log')

    fig10_comps.n(4) = nexttile;
    imagesc((0:delta_t:delta_t*n_fft)*10^15, order, p_injection_cwt);
    colormap turbo;
    hold on 
    yline(harmonic, 'w-')
    xline(tau_pump(delay_idx)*10^15, 'w-')
    xline(tau_probe(delay_idx)*10^15, 'w:')
    caxis([23, 30]);
    title(['$$ \partial_t j_{Injection}, \tau_{Delay}=$$', num2str(delay_between_pulses(delay_idx)), ' fs'])
    xlim([0, 1000])
    xlabel('t in fs')
    ax = gca;
    set(ax, 'YDir', 'normal');
    set(ax, 'YScale', 'log')
    
    cbh = colorbar;
    cbh.Ticks = linspace(23, 30, 10);
    cbh.Layout.Tile = 'east'; 
    STANDARDIZE_FIGURE(fig10_comps);
end


if hue_plot
    figure(4)
    overall_power = log(power_density_harmonic(:, 1:idx+100)) .* exp(1i .* phase_density_harmonic(:, 1:idx+100));
    z_black_to_white_overall = mat2rgbCplx(overall_power, max(max(abs(overall_power)), 1));
    imagesc(abs(overall_power), 'CData', z_black_to_white_overall)
    xlim([idx - 60, idx + 60]);
    title('overall')
    ax = gca;
    ax.YDir = 'normal';
end