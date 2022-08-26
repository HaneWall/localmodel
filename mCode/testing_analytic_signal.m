% testing analyticsignal and phase evaluation 
clear all; 
fc = 20;
delta_t = 0.00001;
t = -1:delta_t:2;
n_fft = length(t);
f = 1/(n_fft*delta_t)*(-(floor((n_fft)/2)):1:(floor((n_fft)/2)));

t_peak = 0.7;
t_ref = 0.7;
width = 0.01;
signal = cos(2*pi*fc*(t-t_peak)) .* exp(-((t-t_peak).^2)/width);
%signal = exp(-((t-t_peak).^2)/width);
[z, phase] = analytic_signal_with_t(signal, t, t_ref);
      
fig1_comps.fig = gcf;
fig1_comps.t1 = tiledlayout(fig1_comps.fig, 1, 2);
fig1_comps.n(1) = nexttile;
p1 = plot(t, signal);
hold on 
p2 = plot(t, abs(z));
set(p2, 'color', 'r', 'LineStyle','-.', 'LineWidth', 1.5);
legend([p1, p2], "Signal", "Hilbert-Envelope")
xlabel("$$t$$ in s")
ylabel("amplitude")
title("$$t_{peak} = $$" + t_peak + "s");

fig1_comps.n(2) = nexttile;
p1 = plot(f, phase);
hold on
p2 = plot(f, fftshift(abs(fft(signal))/max(abs(fft(signal)))));
set(p1, 'color', 'k', 'LineWidth', 1.5);
set(p2, 'color', 'r', 'LineStyle','-.', 'LineWidth', 1.5);
xlim([0, 30])
ylim([-pi, pi])
xlabel("$$f$$ in Hz")
ylabel("$$\angle(\mathcal{F}(signal)(f))$$")
title("$$t_{ref} = $$" + t_ref + "s");
STANDARDIZE_FIGURE(fig1_comps);
