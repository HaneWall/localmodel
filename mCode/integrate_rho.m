%integrate rho_sfi
amplitude = 10^19;
wavelength = 800*10^(-9);
adk_0 = 1e15;
tspan = [0 2*10^(-14)];
rho0 = 0;
[t,rho_sfi] = ode45(@(t,rho_sfi) (1-rho_sfi)*ADK_rate(adk_0, amplitude, wavelength, t), tspan, rho0);
plot(t,rho_sfi,'-o')

