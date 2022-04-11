% parameters to compare
%clear all;
q = 1.60217662*10^(-19);
adk0 = 8.540703969149006e+12;
max_efield = 9508694720.09953;
ionization_potential = 7.5 * q;

e_range = logspace(9, 11, 200);
e_range_perfect = linspace(10^9, max_efield, 200);
gamma = Gamma_ADK(ionization_potential, e_range, 1, 0, 0);
gamma_sim = adk0*(e_range./max_efield).^(13);
gamma_perfect = tangent_Gamma_ADK(e_range_perfect, 7.5);
figure()
loglog(e_range, gamma)
hold on 
loglog(e_range, gamma_sim)
loglog(e_range_perfect, gamma_perfect)
scatter(max_efield, adk0)
xlabel('log{E}');
ylabel('\Gamma_0 in s^{-1}');

