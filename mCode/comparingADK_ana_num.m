% parameters to compare
clear all;
t_end = 1000e-15;
delta_t = 5e-18;
t = 0:delta_t:t_end;
q = 1.60217662*10^(-19);
e_field_pump = gaussian_efield_new(intensity2amplitude(15e16), 2100e-9, 140e-15, 430e-15, t);
max_efield = max(abs(e_field_pump));
ionization_potential = 7.5 * q;

e_range = logspace(9, 11, 200);
e_range_perfect = linspace(10^9, max_efield, 200);
gamma = Gamma_ADK(ionization_potential, e_range, 1, 0, 0);
gamma_perfect = tangent_Gamma_ADK(e_range_perfect, 7.5);

figure(1)
fig1_comps.fig = gcf;
hold on
p1 = plot(log10(e_range), log10(gamma));
p2 = plot(log10(e_range_perfect), log10(gamma_perfect));
p3 = plot(log10(max_efield), log10(max(gamma_perfect)));
hold off

xlabel('$$\log_{10}{E}$$');
ylabel('$$\log_{10}{\Gamma}$$');

set(p1, 'LineStyle', '-', 'LineWidth', 3, 'Color', 'black');
set(p2, 'LineStyle', '--', 'LineWidth', 3, 'Color', 'Red');
set(p3, 'Marker', 'o', 'Color', 'Red', 'Linewidth', 1, 'MarkerSize', 7);
grid on
legend([p1, p2], '$$\Gamma_{ADK}$$', '$$\Gamma$$', 'location', 'southeast')

STANDARDIZE_FIGURE(fig1_comps);
%SAVE_MY_FIGURE(fig1_comps, 'Figures/ComparingGamma.png', 'small');


