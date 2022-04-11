function [gamma] = tangent_Gamma_ADK(efield, bandgap_in_eV)
h = 1e-3;
q = -1.60217662e-19;
e_hat = max(abs(efield));
gamma_hat = Gamma_ADK(bandgap_in_eV*(-q), e_hat, 1, 0, 0);
a = e_hat./gamma_hat .* (Gamma_ADK(bandgap_in_eV*(-q), e_hat+h, 1, 0, 0) - Gamma_ADK(bandgap_in_eV*(-q), e_hat, 1, 0, 0))./h ;
gamma = gamma_hat .* abs(efield./e_hat).^a ;
end
