%permutation plot that shows dampeneing for large omega behaves like multinomial distribution
clear all;
c = 299792458;
q = 1.60217662e-19;
me = 9.10938e-31;
n0 = 2.2e28; 
bandgaps = [7.7];
wavelength_probe = 800e-9;
wavelength_pump = 2100e-9;


i_pump  = 11.5e16;
i_probe = 0.005e16;

amplitude_pump = intensity2amplitude(i_pump);
amplitude_probe = intensity2amplitude(i_probe);

% integration params
t_end = 1000e-15;
delta_t = 5e-18;
t = 0:delta_t:t_end;
L = length(t);

%delay times
tau_pump = 500e-15;
tau_probe = 500e-15;

e_field_probe = zeros(3, L);
e_field_probe(1,:) = gaussian_efield_new(amplitude_probe, wavelength_probe, 75e-15, tau_probe, t);
normed_e_field = zeros(1, L);

third_term_parall = zeros(3, L);
third_term_perpe = zeros(3, L);

e_field_pump_parall = zeros(3, L);
e_field_pump_perpe = zeros(3, L);

injection = zeros(2, L);
brunel = zeros(2, L);
kerr = zeros(2, L);


n_fft = 2^nextpow2(L); %zero padding
f = 1/delta_t*(0:(n_fft/2))/n_fft;
delta_f = 1/delta_t;
f_pump = c / 2.1e-06;
f_probe = c / 8e-07;
first_harm = 2*f_pump + f_probe;
[~, idx] = min(abs(first_harm - f));
[~, idx_pump] = min(abs(f_pump - f));
[~, idx_probe] = min(abs(f_probe - f));

for i =1:2
    if i == 1
        % parall case 
        amplitude_sum = amplitude_pump + amplitude_probe;
        e_field_pump_parall(1,:) = gaussian_efield_new(amplitude_pump, wavelength_pump, 140e-15, tau_pump, t);
        e_field = e_field_pump_parall + e_field_probe;
        for j = 1:L
            normed_e_field(:,j) = norm(e_field(:,j));
        end
        displacements_x = displacement_x_new(bandgaps, normed_e_field + 1, e_field);
        ADK = tangent_Gamma_ADK(normed_e_field, bandgaps);
        %rho_sfi = integrate_population_cb(ADK, delta_t, t);
        rho_sfi = integrate_cb_wo_sat(ADK, delta_t, t);
        drho = gradient(rho_sfi, delta_t);
        third_term_parall(1, :) = gradient(displacements_x(1,:).*drho, delta_t);
        v0 = 0;
        brunel_vec = n0 * q * q/me * e_field.*rho_sfi; 
        kerr_vec = kerr_current_binomial(e_field_pump_parall, e_field_probe, delta_t);
        injection_vec =  n0 .* q .* third_term_parall; 
        
        brunel(i, :) = brunel_vec(1, :);
        kerr(i, :) = kerr_vec(1, :);
        injection(i, :) = injection_vec(1, :);
    else
        % perpendicular case
        e_field_pump_perpe(2,:) = gaussian_efield_new(amplitude_pump, wavelength_pump, 140e-15, tau_pump, t);
        e_field = e_field_pump_perpe + e_field_probe;
        for j = 1:L
            normed_e_field(:,j) = norm(e_field(:,j));
        end
        displacements_x = displacement_x_new(bandgaps, normed_e_field + 1, e_field);
        ADK = tangent_Gamma_ADK(normed_e_field, bandgaps);
        %rho_sfi = integrate_population_cb(ADK, delta_t, t);
        rho_sfi = integrate_cb_wo_sat(ADK, delta_t, t);
        drho = gradient(rho_sfi, delta_t);
        third_term_perpe(1,:) = gradient(displacements_x(1,:).*drho, delta_t);
        v0 = 0;
        brunel_vec_perpe = n0 * q * q/me * e_field.*rho_sfi;
        kerr_vec_perpe = kerr_current_binomial(e_field_pump_parall, e_field_probe, delta_t);
        injection_vec_perpe =  n0 .* q .* third_term_perpe;
        
        brunel(i, :) = brunel_vec_perpe(1, :);
        kerr(i, :) = kerr_vec_perpe(1, :);
        injection(i, :) = injection_vec_perpe(1, :);
        
    end

end

ft_brunel_para = fft(brunel(1,:), n_fft);
ft_kerr_para = fft(kerr(1,:), n_fft);
ft_injection_para = fft(injection(1,:), n_fft);

ft_brunel_perp = fft(brunel(2,:), n_fft);
ft_kerr_perp = fft(kerr(2,:), n_fft);
ft_injection_perp = fft(injection(2,:), n_fft);

spec_brunel_para = abs(ft_brunel_para(1:n_fft/2+1)/n_fft);
spec_kerr_para = abs(ft_kerr_para(1:n_fft/2+1)/n_fft);
spec_injection_para = abs(ft_injection_para(1:n_fft/2+1)/n_fft);

spec_brunel_perp = abs(ft_brunel_perp(1:n_fft/2+1)/n_fft);
spec_kerr_perp = abs(ft_kerr_perp(1:n_fft/2+1)/n_fft);
spec_injection_perp = abs(ft_injection_perp(1:n_fft/2+1)/n_fft);

slope = spec_injection_para(idx)./spec_injection_perp(idx);
%m = round(spec_injection_para(idx)./spec_injection_perp(idx));
m_injection = round(slope);
m_brunel = round(slope) + 2;
n_injection = 1:1:((m_injection-1)/2);
n_brunel = 1:1:((m_brunel-1)/2);

n_pu_i = 2*n_injection;
n_pu_b = 2*n_brunel;

n_pr = 1;
permutations_i = zeros(length(n_injection),1);
permutations_b = zeros(length(n_brunel),1);

for j=1:length(n_injection)
    permutations_i(j) = multinomial_degen(m_injection, n_pu_i(j), n_pr) * (2*n_injection(j)*idx_pump + idx_probe);
end

for j=1:length(n_brunel)
    permutations_b(j) = multinomial_degen(m_brunel, n_pu_b(j), n_pr) * (2*n_brunel(j)*idx_pump + idx_probe)^(-1);
end

figure(1)
fig1_comps.fig = gcf;
fig1_comps.t1 = tiledlayout(fig1_comps.fig, 2, 1);
fig1_comps.n(1) = nexttile;
p4 = plot(f./f_pump, spec_injection_para./spec_injection_para(idx));
hold on
p1 = stem(n_pu_i + f_probe/f_pump, permutations_i/max(permutations_i));
p2 = stem(n_pu_i + f_probe/f_pump, permutations_i/(max(permutations_i)*m_injection));
p3 = plot(f./f_pump, spec_injection_perp./spec_injection_para(idx));
set(p1, 'LineWidth', 1, 'color', 'blue');
set(p2, 'Linewidth', 1, 'color', 'red');
set(p3, 'LineWidth', 1, 'color', 'red');
set(p4, 'Linewidth', 1, 'color', 'blue');
xlim([0, 20]);
ylim([10e-11, 10e2]);
a = gca;
set(a, 'YScale', 'log');
title('Injection')
legend([p4, p3], '$$\mathcal{F}(\partial_t \mathbf{j}_{\|})$$', '$$\mathcal{F}(\partial_t \mathbf{j}_{\perp})$$', 'location', 'northeast')
%STANDARDIZE_FIGURE(fig_comps1)

peaks = [];
for k = 2:2:m_injection-1
    peaks = [peaks, spec_injection_para(k*idx_pump + idx_probe)];
end
peaks_normed = peaks./max(peaks);
per_normed = permutations_i./max(permutations_i);

ratio = peaks_normed./per_normed';


fig1_comps.n(2) = nexttile;
p4 = plot(f./f_pump, spec_brunel_para./spec_brunel_para(idx));
hold on
p1 = stem(n_pu_b + f_probe/f_pump, permutations_b/max(permutations_b));
p2 = stem(n_pu_b + f_probe/f_pump, permutations_b/(max(permutations_b)*m_brunel));
p3 = plot(f./f_pump, spec_brunel_perp./spec_brunel_para(idx));
set(p1, 'LineWidth', 1, 'color', 'blue');
set(p2, 'Linewidth', 1, 'color', 'red');
set(p3, 'LineWidth', 1, 'color', 'red');
set(p4, 'Linewidth', 1, 'color', 'blue');
xlim([0, 20]);
ylim([10e-11, 10e2]);
fig1_comps.tile2.plotXLabel = xlabel('$$\omega / \omega_{pump}$$');
a = gca;
set(a, 'YScale', 'log');
title('Brunel')
legend([p4, p3], '$$\mathcal{F}(\partial_t \mathbf{j}_{\|})$$', '$$\mathcal{F}(\partial_t \mathbf{j}_{\perp})$$', 'location', 'northeast')
STANDARDIZE_FIGURE(fig1_comps)





