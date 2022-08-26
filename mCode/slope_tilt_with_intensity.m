%omega_probe in respect to the delay (parallel polarization)
clear all;
%physical constants
c = 299792458;
q = 1.60217662e-19;
me = 9.10938e-31;
n0 = 2.2e28;                                                 %molecular density for si02

%plot parameters
linear_regression = true;
plot_overall_over_time = false;

%simulation parameters 
bandgap = 7.5; % in eV

amplitude_probe= intensity2amplitude(0.015e16);  %0.015 TWcm^-2
%amplitude_pump = intensity2amplitude(12e16);     %12 TWcm^-2
I_length = 5;
intensity_range = linspace(5e16, 20e16, I_length);
amplitude_pump = intensity2amplitude(intensity_range);     %13 TWcm^-2

wavelength_probe = 800e-9;
wavelength_pump = 2100e-9;

%integration parameters
t_end = 1000e-15;
delta_t = 3.8e-18;
t = 0:delta_t:t_end;

%allocate memory for first harmonic
L = length(t);
n_fft = 2^nextpow2(L); %zero padding

t_all = delta_t:delta_t:n_fft*delta_t;
t_all_fs = t_all*10^15;

f = 1/delta_t*(0:(n_fft/2))/n_fft;
delta_f = 1/delta_t;
f_pump = c / 2.1e-06;
f_probe = c / 8e-07;
first_harm = 2*f_pump + f_probe;
[~, idx] = min(abs(first_harm - f));
[~, idx_pump] = min(abs(f_pump - f));
[~, idx_probe] = min(abs(f_probe - f));

%rectangle window
width_rectangle_freq = 80;

%shifted parameters, need two indices, because i have negative and positive
%frequencies 
f_shifted = (-n_fft/2:n_fft/2 - 1) .* delta_f/n_fft;
middle = int64(n_fft/2);
[~, first_idx_s] = min(abs(first_harm + f_shifted(1:middle)));
%[~, second_idx_s_proto] = min(abs(first_harm + f_probe - f_shifted(middle+1:end)));
second_idx_s = first_idx_s + 2*(middle-first_idx_s);
%second_idx_s = middle + second_idx_s_proto;


%define window function for ectracting only the first harmonic 
%rectangle version
rectangle = zeros(1, n_fft);
rectangle(first_idx_s-width_rectangle_freq/2:first_idx_s+width_rectangle_freq/2) = 1;
%rectangle(second_idx_s-width_rectangle_freq/2:second_idx_s+width_rectangle_freq/2) = 1;

% FWHM in terms of intensity 
fwhm_pump = 140e-15;
fwhm_probe = 45e-15;


%delay times
delay_between_pulses = -150e-15:2.5e-15:150e-15;
N = length(delay_between_pulses);
tau_pump = 500e-15;
tau_probe = tau_pump + delay_between_pulses;

kerr_container = zeros(length(delay_between_pulses), n_fft);
injection_container = zeros(length(delay_between_pulses), n_fft);
brunel_container = zeros(length(delay_between_pulses), n_fft);
overall_container = zeros(I_length, length(delay_between_pulses), n_fft);

%center of mass evaluation of t-axes --> maybe linear dependence 
A = zeros(I_length, N, 1);
center = zeros(I_length, N, 1);
center_of_mass_t = zeros(I_length, N, 1);

max_idx_injection = zeros(I_length, N, 1);
max_idx_kerr = zeros(I_length, N, 1);
max_idx_brunel = zeros(I_length, N, 1);
max_idx_overall = zeros(I_length, N, 1);

for I_idx = 1:1:I_length
    for N_idx = 1:1:N
        e_field_probe = zeros(3, L);
        e_field_pump = zeros(3, L);
        third_term = zeros(3, L);
        normed_e_field = zeros(1, L);
        
        %run simulation 
        e_field_pump(1, :) = gaussian_efield_new(amplitude_pump(I_idx), wavelength_pump, fwhm_pump, tau_pump, t);
        e_field_probe(1, :) = gaussian_efield_new(amplitude_probe, wavelength_probe, fwhm_probe, tau_probe(N_idx), t);
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
            
        ft_plasma_current = fft(plasma_current_density(1, :), n_fft);
        ft_brunel_current = fft(brunel_current_density(1, :), n_fft);
        ft_injection_current = fft(injection_current_density(1, :), n_fft);
        ft_kerr_current = fft(kerr(1, :), n_fft); 
        ft_overall_current = fft(overall_current_density(1, :), n_fft);
        
        real_ft_brunel =fftshift(ft_brunel_current);
        real_ft_injection = fftshift(ft_injection_current);
        real_ft_kerr = fftshift(ft_kerr_current);
        real_ft_overall = fftshift(ft_overall_current);
        
        extracted_harm_brunel = real_ft_brunel.*rectangle;
        extracted_harm_injection = real_ft_injection.*rectangle;
        extracted_harm_kerr = real_ft_kerr.*rectangle;
        extarcted_harm_overall = real_ft_overall.*rectangle;
        
        inverse_brunel = ifft(extracted_harm_brunel);
        inverse_injection = ifft(extracted_harm_injection);
        inverse_kerr = ifft(extracted_harm_kerr);
        inverse_overall = ifft(extarcted_harm_overall);
        
        overall_container(I_idx, N_idx, :) = abs(inverse_overall);
        
        [~, max_idx_injection(I_idx, N_idx)] = max(abs(inverse_injection));
        [~, max_idx_brunel(I_idx, N_idx)] = max(abs(inverse_brunel));
        [~, max_idx_kerr(I_idx, N_idx)] = max(abs(inverse_kerr));
        [~, max_idx_overall(I_idx, N_idx)] = max(abs(inverse_overall));
        
        A(I_idx, N_idx) = trapz((t_all_fs(1:ceil(n_fft/2))), overall_container(I_idx, N_idx, 1:ceil(n_fft/2)));
        center(I_idx, N_idx) = trapz(t_all_fs(1:ceil(n_fft/2)), t_all_fs(1:ceil(n_fft/2)) .* squeeze(overall_container(I_idx, N_idx, 1:ceil(n_fft/2)))');
    end
end

center_of_mass_t = center./A;

if linear_regression
    slope_inj=zeros(I_length, 1);
    slope_brunel=zeros(I_length, 1);
    slope_kerr=zeros(I_length, 1);
    slope_overall=zeros(I_length, 1);
    X = [ones(N,1) transpose(delay_between_pulses)*10^15];
    for I_idx = 1:1:I_length
        b_inj= X\(max_idx_injection(I_idx, :)'*delta_t*10^15);
        slope_inj(I_idx) = b_inj(2);

        b_brunel = X\(max_idx_brunel(I_idx, :)'*delta_t*10^15);
        slope_brunel(I_idx) = b_brunel(2);

        b_kerr = X\(max_idx_kerr(I_idx, :)'*delta_t*10^15);
        slope_kerr(I_idx) = b_kerr(2);

        b_overall = X\(max_idx_overall(I_idx, :)'*delta_t*10^15);
        slope_overall(I_idx) = b_overall(2);
    end

    figure(2)
    fig2_comps.fig = gcf;
    p1 = scatter(intensity_range*10^-16, slope_kerr);
    hold on
    p2 = scatter(intensity_range*10^-16, slope_inj);
    p3 = scatter(intensity_range*10^-16, slope_brunel);
    p4 = scatter(intensity_range*10^-16, slope_overall);
    
    ylabel("tilt")
    xlabel("$$I$$ in $$TWcm^(-2)$$")
    legend([p1, p2, p3, p4], "Kerrtilt", "Injectiontilt", "Bruneltilt", "Overalltilt")
    STANDARDIZE_FIGURE(fig2_comps);
end                     
    

if plot_overall_over_time
    figure(10)
    fig10_comps.fig = gcf;
    c1 = imagesc(t_all_fs, delay_between_pulses*10^15, log10(overall_container));
    colormap jet;
    ax = gca;
    xlim([300, 700]);
    caxis([25.7, 27.7]);
    set(ax, 'YDir', 'normal');
    hold on 
    p1=plot(t_all_fs, (t_all_fs-500) .* (1/b_kerr(2)));
    p2=plot(t_all_fs, (t_all_fs-500) .* (1/b_inj(2)));
    p3=scatter(max_idx_overall .* delta_t * 10^15, delay_between_pulses*10^15);
    p4=scatter(center_of_mass_t, delay_between_pulses*10^15);
    set(p1, 'LineWidth', 1.5, 'LineStyle', '-.', 'Color', 'w');
    set(p2, 'LineWidth', 1.5, 'LineStyle', '-', 'Color', 'w');
    set(p3, 'LineWidth', 1, 'MarkerEdgeColor', 'w');
    set(p4, 'LineWidth', 1, 'MarkerEdgeColor', 'w', 'Marker','square');
    title("$$\log_{10}\|\mathcal{F}^{-1}(\partial_t j(\omega_1))\|$$")
    legend([p1, p2, p3, p4], "Kerrtilt", "Injectiontilt", "PeakMaxima $$\omega_1$$", "center of mass $$\omega_1$$")
    ylabel("$$\tau_{delay}$$ in fs")
    xlabel("$$t$$ in fs")
    STANDARDIZE_FIGURE(fig10_comps);
end
