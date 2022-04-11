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
figure()
loglog(e_range, gamma)
hold on 
loglog(e_range_perfect, gamma_perfect)
scatter(max_efield, max(gamma_perfect))
xlabel('log{E}');
ylabel('\Gamma_0 in s^{-1}');

