%script that does compare the first harmonic n = 1 for 2*n omega_pump +
%omega_probe in respect to the delay (parallel polarization)
clear all;
%physical constants
c = 299792458;
q = 1.60217662e-19;
me = 9.10938e-31;
n0 = 2.2e28;                                     %molecular density for si02

phase_plot=true;
hue_plot=false;
delay_plot_all=false;
delay_filter_plot=false;

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
delay_between_pulses = linspace(-200e-15, 200e-15, 100); %80
tau_pump = ones(1, 100) * 500e-15;
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

for i = 1:length(delay_between_pulses)
    e_field_pump = gaussian_efield_new(amplitude_pump, wavelength_pump, fwhm_pump, tau_pump(i), t);
    e_field_probe = gaussian_efield_new(amplitude_probe, wavelength_probe, fwhm_probe, tau_probe(i), t);
    e_field = e_field_pump + e_field_probe;
    displacements_x = displacement_x_new(bandgap, abs(e_field) + 1, e_field);
    ADK = tangent_Gamma_ADK(e_field, bandgap);
    
    rho_sfi = integrate_population_cb(ADK, delta_t, t);
    drho = gradient(rho_sfi, delta_t);
    third_term = gradient(displacements_x.*drho, delta_t);
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
    idx_half_band = 55;
    modified_log_spec = log_power_spec;
    for i=0:7
        modified_log_spec(:, (2*i+1)*idx_pump - idx_half_band:(2*i+1)*idx_pump + idx_half_band) = 30; 
        modified_log_spec(:, (2*i+1)*idx_probe - idx_half_band:(2*i+1)*idx_probe + idx_half_band) = 30;
    end
    figure(6)
    % Get the handle of figure(n).
    fig1_comps.fig = gcf;
    p1 = imagesc(f./f_pump, delay_between_pulses.*10^(15), modified_log_spec, [47.5, 52]);
    fig1_comps.tile1.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig1_comps.tile1.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), modified_log_spec, [49 50 51]);
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
    figure(7)
    % Get the handle of figure(n).
    fig2_comps.fig = gcf;
    imagesc(f./f_pump, delay_between_pulses.*10^(15), log_power_spec, [47.5, 52]);
    fig2_comps.tile1.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig2_comps.tile1.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
    title('Overall')
    xline(f(idx)/f_pump, 'w-');
    xline(f(idx_pump)/f_pump, 'w-.');
    xline(f(idx_probe)/f_pump, 'w--');
    xlim([0, 11]);
    ax = gca;
    ax.YDir = 'normal';
    colormap jet;
    STANDARDIZE_FIGURE(fig2_comps);
    
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
    title('Overall')
    xline(f(idx)/f_pump, 'w-');
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    ax = gca;
    ax.YDir = 'normal';
    fig3_comps.tile1.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig3_comps.tile1.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');

    

    %STANDARDIZE_FIGURE(fig3_comps);

    fig3_comps.n(2) = nexttile;
    fig3_comps.tile2.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig3_comps.tile2.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
    imagesc(f/f_pump, delay_between_pulses.*10^(15), log_kerr_spec, [47.5, 52]);
    colormap jet;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_kerr_spec, [49 50 51]);
    c1.LineColor = 'black';
    c1.LineWidth = 0.5;
    caxis(cRange)
    title('Kerr')
    xline(f(idx)/f_pump, 'w-');
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    ax = gca;
    ax.YDir = 'normal';
    fig3_comps.tile2.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig3_comps.tile2.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
    
    fig3_comps.n(3) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), log_injection_spec, [47.5, 52]);
    colormap jet;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_injection_spec, [49 50 51]);
    c1.LineColor = 'black';
    c1.LineWidth = 0.5;
    caxis(cRange)
    title('Injection')
    xline(f(idx)/f_pump, 'w-');
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    ax = gca;
    ax.YDir = 'normal';
    fig3_comps.tile3.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig3_comps.tile3.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
    
    
    fig3_comps.n(4) = nexttile;
    imagesc(f/f_pump, delay_between_pulses.*10^(15), log_brunel_spec, [47.5, 52]);
    colormap jet;
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f/f_pump, delay_between_pulses.*10^(15), log_brunel_spec, [49 50 51]);
    c1.LineColor = 'black';
    c1.LineWidth = 0.5;
    colormap jet;
    caxis(cRange)
    title('Brunel')
    xline(f(idx)/f_pump, 'w-');
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    ax = gca;
    ax.YDir = 'normal';
    fig3_comps.tile4.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
    fig3_comps.tile4.plotYLabel = ylabel('$$\tau_{Delay}$$ in fs');
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
    title('Phase zw. Injection & Brunel')
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
    caxis(cRange)
    xlim([f(idx - 60)/f_pump, f(idx + 60)/f_pump]);
    xline(f(idx)/f_pump, 'b-');
    title('Phase zw. Injection & Kerr')
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
    title('Phase zw. Brunel & Kerr')
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