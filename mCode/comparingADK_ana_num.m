% parameters to compare
%clear all;
q = 1.60217662*10^(-19);
adk0 = 8.540703969149006e+12;
%adk0 = 8.540703969149006e+12 / 2;
max_efield = 9508694720.09953;
ionization_potential = 7.5 * q;

e_range = logspace(9, 11, 200);
gamma = Gamma_ADK(ionization_potential, e_range, 1, 0, 0);
gamma_sim = adk0*(e_range./max_efield).^(13);
figure()
loglog(e_range, gamma)
hold on 
loglog(e_range, gamma_sim)
scatter(max_efield, adk0)
xlabel('log{E}');
ylabel('\Gamma_0 in s^{-1}');

