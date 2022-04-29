clear all;
%physical constants
q = -1.60217662e-19;
me = 0.635 * 9.10938e-31;
n0 = 2.2e28; %molecular density for si02
c = 299792458;
bandgaps = [7.5];

%simualtion parameters 
%amplitude_probe= intensity2amplitude(0.015e16); %0.015 TWcm^-2
amplitude_probe = 0; %0.015 TWcm^-2
amplitude_pump = intensity2amplitude(12e16); %12 TWcm^-2
amplitude_sum = amplitude_pump + amplitude_probe;
wavelength_probe = 800e-9;
wavelength_pump = 2100e-9;

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
ADK = tangent_Gamma_ADK(e_field, bandgaps);

rho_sfi = integrate_population_cb(ADK, delta_t, t);

%compute time derivative plasma current density 
drho = gradient(rho_sfi, delta_t);
third_term = gradient(displacements_x.*drho, delta_t);
v0 = 0;
plasma_current_density = n0 * q * (q/me * e_field.*rho_sfi + v0*drho+ third_term);
brunel_current_density = n0 * q * q/me * e_field.*rho_sfi;
kerr = kerr_current_density(e_field, delta_t); 
injection_current_density =  n0 * q * third_term; 
overall_current_density = plasma_current_density + kerr ;

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

figure(1);
% Get the handle of figure(n).
fig1_comps.fig = gcf;
fig1_comps.t1 = tiledlayout(fig1_comps.fig, 1, 2);
fig1_comps.n(1) = nexttile;
hold on
p1 = area(t - tau_pump, ADK./max(ADK).*max(rho_sfi));
p2 = plot(t - tau_pump, rho_sfi, '-');
xlim([-1e-13, 1e-13])
fig1_comps.tile1.plotXLabel = xlabel('$$t-\tau_{delay}$$ in s');
legend([p1, p2], '$$\Gamma$$', '$$\rho$$');
set(p2, 'LineWidth', 1, 'Color', 'Blue');
set(p1, 'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
%STANDARDIZE_FIGURE(fig1_comps);

%figure(2);
%fig2_comps.fig = gcf;
fig1_comps.n(2) = nexttile;
hold on
p3 = plot(omegas/omega_pump,log10(P_brunel_current(1:n_fft/2 + 1)/max(P_overall_current)));
p4 = plot(omegas/omega_pump,log10(P_injection_current(1:n_fft/2 + 1)/max(P_overall_current)));
p5 = plot(omegas/omega_pump,log10(P_kerr_current(1:n_fft/2 + 1)/max(P_overall_current)));
p6 = plot(omegas/omega_pump,log10(P_overall_current(1:n_fft/2 + 1)/max(P_overall_current)));
for k = 1:20
    %xline((2*k*omega_pump + omega_probe),'-');
    xline(k, '-.');
    %xline(k*omega_probe, ':');
end
set(p3, 'LineWidth', 1);
set(p4, 'LineWidth', 1);
set(p5, 'LineWidth', 1);
set(p6, 'LineWidth', 1);
xlim([0, 12])
legend([p3, p4, p5, p6],'Brunel','Injection', 'Kerr', 'Overall', 'location', 'southeast');
fig1_comps.tile2.plotxLabel = xlabel('Harmonische');
fig1_comps.tile2.plotyLabel = ylabel('$$\log_{10}|\mathcal{F}(\partial_t j)|^2$$');
STANDARDIZE_FIGURE(fig1_comps);



