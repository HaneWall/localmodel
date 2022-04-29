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

figure(1);
% Get the handle of figure(n).
fig1_comps.fig = gcf;
hold on 
m = sqrt(parall_overall./perpe_overall);
p1 = plot(i_pump_ranges, m);
p2 = plot(i_pump_ranges, sqrt(parall_kerr./perpe_kerr));
p3 = plot(i_pump_ranges, sqrt((parall_kerr + parall_brunel)./(perpe_kerr + perpe_brunel)));
fig1_comps.tile1.plotXLabel = xlabel('$$I_{pump}$$ in $$Wm^{-2}$$');
fig1_comps.tile1.plotYLabel = ylabel('$$m = \sqrt{I_{\mid\mid}/I_{\perp}}$$');
%title('$$E_{Gap}=7.5 eV , \Delta t = 1e-18s, \tau_{delay} = 0$$', 'Interpreter','latex')
set(p1, 'LineWidth', 1, 'Marker','o');
set(p2, 'LineWidth', 1);
set(p3, 'LineWidth', 1, 'Marker','o');
legend('m', 'Kerr', 'Kerr+Brunel', 'location', 'northwest')
ylim([2,14])
grid on
xlim([i_pump_ranges(1), i_pump_ranges(end)])
STANDARDIZE_FIGURE(fig1_comps)
