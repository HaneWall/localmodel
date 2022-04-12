% creating m figure
load("parall_kerr.mat", "parall_kerr");
load("parall_brunel.mat", "parall_brunel");
load("parall_injection.mat", "parall_injection");
load("parall_overall.mat", "parall_overall");

load("perpe_kerr.mat", "perpe_kerr");
load("perpe_brunel.mat", "perpe_brunel");
load("perpe_injection.mat", "perpe_injection");
load("perpe_overall.mat", "perpe_overall");

%load("intensity_ranges.mat", "e_pump_ranges")
i_pump_ranges = linspace(3e16, 18e16, 20);

m = sqrt(parall_overall./perpe_overall);
plot(i_pump_ranges, m)
hold on 
plot(i_pump_ranges, sqrt(parall_kerr./perpe_kerr))
plot(i_pump_ranges, sqrt((parall_kerr + parall_brunel)./(perpe_kerr + perpe_brunel)))
xlabel('I_{pump} in W/m^2', 'interpreter', 'tex')
ylabel('$$m = \sqrt{I_{\mid\mid}/I_{\perp}}$$', 'Interpreter','latex')
title('$$E_{Gap}=7.5 eV , \Delta t = 1e-18s, \tau_{delay} = 0$$', 'Interpreter','latex')
legend('m', 'Kerr', 'Kerr+Brunel', 'Kerr+Brunel+Injection')
