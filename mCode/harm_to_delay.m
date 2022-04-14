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
delay_plot_all=true;

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
log_power_spec = log(power_density_harmonic);
log_kerr_spec = log(power_kerr_harmonic);
log_injection_spec = log(power_injection_harmonic);
log_brunel_spec = log(power_brunel_harmonic);
f = 1/delta_t*(0:(n_fft/2))/n_fft;
delta_f = 1/delta_t;
f_pump = c / 2.1e-06;
f_probe = c / 8e-07;
first_harm = 2*f_pump + f_probe;
[~, idx] = min(abs(first_harm - f));
[~, idx_pump] = min(abs(f_pump - f));
[~, idx_probe] = min(abs(f_probe - f));


if delay_plot_all
    subplot(3,4,[1, 2, 3, 4]);
    imagesc(f, delay_between_pulses, log_power_spec, [111, 120]);
    hold on
    title('Overall')
    xline(f(idx), 'w-');
    xline(f(idx_pump), 'w-.');
    xline(f(idx_probe), 'w--');
    xlim([0,f(idx+1500)]);
    ax = gca;
    ax.YDir = 'normal';
    
    subplot(3,4,5);
    imagesc(f, delay_between_pulses, log_power_spec, [111, 120]);
    hold on
    title('Overall')
    xline(f(idx), 'w-');
    xline(f(idx_pump), 'w-.');
    xline(f(idx_probe), 'w--');
    xlim([0,f(idx+1500)]);
    ax = gca;
    ax.YDir = 'normal';

    subplot(3,4,6);
    imagesc(f, delay_between_pulses, log_kerr_spec, [111, 120]);
    hold on
    title('Kerr')
    xline(f(idx), 'w-');
    xline(f(idx_pump), 'w-.');
    xline(f(idx_probe), 'w--');
    xlim([0,f(idx+1500)]);
    ax = gca;
    ax.YDir = 'normal';
    
    subplot(3,4,7);
    imagesc(f, delay_between_pulses, log_injection_spec, [111, 120]);
    hold on
    title('Injection')
    xline(f(idx), 'w-');
    xline(f(idx_pump), 'w-.');
    xline(f(idx_probe), 'w--');
    xlim([0,f(idx+1500)]);
    ax = gca;
    ax.YDir = 'normal';
    
    
    subplot(3,4,8);
    imagesc(f, delay_between_pulses, log_brunel_spec, [111, 120]);
    hold on
    title('Brunel')
    xline(f(idx), 'w-');
    xline(f(idx_pump), 'w-.');
    xline(f(idx_probe), 'w--');
    xlim([0,f(idx+1500)]);
    ax = gca;
    ax.YDir = 'normal';
    
    subplot(3,4,9);
    imagesc(f, delay_between_pulses, phase_density_harmonic);
    colormap jet
    hold on
    title('Overall')
    xline(f(idx), 'w-');
    xline(f(idx_pump), 'w-.');
    xline(f(idx_probe), 'w--');
    xlim([0,f(idx+1500)]);
    ax = gca;
    ax.YDir = 'normal';
    
    subplot(3,4,10);
    imagesc(f, delay_between_pulses, phase_kerr_harmonic);
    colormap jet
    hold on
    title('Kerr')
    xline(f(idx), 'w-');
    xline(f(idx_pump), 'w-.');
    xline(f(idx_probe), 'w--');
    xlim([0,f(idx+1500)]);
    ax = gca;
    ax.YDir = 'normal';
    
    subplot(3,4,11);
    imagesc(f, delay_between_pulses, phase_injection_harmonic);
    colormap jet
    hold on
    title('Injection')
    xline(f(idx), 'w-');
    xline(f(idx_pump), 'w-.');
    xline(f(idx_probe), 'w--');
    xlim([0,f(idx+1500)]);
    ax = gca;
    ax.YDir = 'normal';
    
    subplot(3,4,12);
    imagesc(f, delay_between_pulses, phase_brunel_harmonic);
    colormap jet
    hold on
    title('Brunel')
    xline(f(idx), 'w-');
    xline(f(idx_pump), 'w-.');
    xline(f(idx_probe), 'w--');
    xlim([0,f(idx+1500)]);
    ax = gca;
    ax.YDir = 'normal';
end

if phase_plot
    figure(2)
    subplot(1,3,1)
    imagesc(f, delay_between_pulses, mod(phase_injection_harmonic - phase_brunel_harmonic + pi,2*pi) - pi)
    cRange = caxis; % save the current color range
    hold on 
    [M1,c1] = contour(f, delay_between_pulses, log_power_spec, [111 113 115 117 119]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f(idx - 60), f(idx + 60)]);
    xline(f(idx), 'b-');
    title('phase relation injection and brunel')
    colormap hsv
    ax = gca;
    ax.YDir = 'normal';
    subplot(1,3,2)
    imagesc(f, delay_between_pulses, mod(phase_injection_harmonic - phase_kerr_harmonic + pi,2*pi) - pi)
    hold on 
    [M1,c1] = contour(f, delay_between_pulses, log_power_spec, [111 113 115 117 119]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f(idx - 60), f(idx + 60)]);
    xline(f(idx), 'b-');
    title('phase relation injection and kerr')
    colormap hsv
    ax = gca;
    ax.YDir = 'normal';
    subplot(1,3,3)
    imagesc(f, delay_between_pulses, mod(phase_kerr_harmonic - phase_brunel_harmonic + pi,2*pi) - pi)
    hold on 
    [M1,c1] = contour(f, delay_between_pulses, log_power_spec, [111 113 115 117 119]);
    c1.LineColor = 'white';
    c1.LineWidth = 1;
    caxis(cRange)
    xlim([f(idx - 60), f(idx + 60)]);
    xline(f(idx), 'b-');
    title('phase relation kerr and brunel')
    colormap hsv
    ax = gca;
    ax.YDir = 'normal';
end

if hue_plot
    figure(3)
    overall_power = log(power_density_harmonic(:, 1:idx+100)) .* exp(1i .* phase_density_harmonic(:, 1:idx+100));
    z_black_to_white_overall = mat2rgbCplx(overall_power, max(max(abs(overall_power)), 1));
    imagesc(abs(overall_power), 'CData', z_black_to_white_overall)
    xlim([idx - 60, idx + 60]);
    title('overall')
    ax = gca;
    ax.YDir = 'normal';
end