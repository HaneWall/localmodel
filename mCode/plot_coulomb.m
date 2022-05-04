 figure(1)
 % Get the handle of figure(n).
 fig2_comps.fig = gcf;
 x = -3:0.01:3;
 p1 = plot(x, -1./abs(x) - 6.*x);
 %p2 = plot(x, - 9.*x);
 ylim([-15, 9]);
 xlim([-2,2]);
 hold on;
 p2 = plot(x, -6.*x);
 set(p1, 'LineWidth', 1.5);
 set(p2, 'LineWidth', 1.5);
 fig2_comps.tile1.plotXLabel = xlabel('$$x$$');
 fig2_comps.tile1.plotYLabel = ylabel('Potential');
 STANDARDIZE_FIGURE(fig2_comps);