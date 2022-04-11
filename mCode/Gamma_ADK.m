function Gamma=Gamma_ADK(I_p, E, Z,l,m)
% Input: 
% I_p : ionization potential / Bandgap
% E: electric field [V/m]
% Z: final charge state (default 1)
% l: angular momentum quantum number (default 0)
% m: magnetic quantum number (default 0)

% Output: 
% Gamma: ionization rate [1/sec]

% calculation below in atomic units (au)

e_0 = 1.6021766208e-19;

I_p_au=I_p./e_0./27.21139609; %au in eV
kappa=sqrt(2*I_p_au);
ns=Z./kappa;
ls=ns-1;
F=E./5.14220652e11;    %E field in au

A_lm=((2.*l+1).*factorial(l+abs(m)))./(2.^(abs(m)).*factorial(abs(m)).*factorial(l-abs(m)));

C_nsls_abs_squared=2.^(2.*ns)./(ns.*gamma(ns+ls+1).*gamma(ns-ls));

Gamma_au=C_nsls_abs_squared.*A_lm.*I_p_au...
    .*(2.*kappa^3./F).^(2.*ns-abs(m)-1).*exp(-(2.*kappa^3./(3.*F)));
% Gamma_au=sqrt(6/pi)*C_nsls_abs_squared.*A_lm.*I_p_au...
%     .*(2.*kappa^3./F).^(2.*ns-abs(m)-3/2).*exp(-(2.*kappa^3./(3.*F)));

Gamma=Gamma_au./24.18884336e-18;
Gamma(isnan(Gamma))=0;
Gamma(isinf(Gamma))=0;
end