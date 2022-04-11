function [der_j_kerr] = kerr_current_density_vec(E_1, E_2, delta_t)
    %UNTITLED7 Summary of this function goes here
    %   Detailed explanation goes here
    eps0 = 8.8541878128e-12;
    % third order chi SiO2
    chi_3 = 1.94e-22;
    L = length(E_1);
    p = zeros(3, L);
    current_dens_vec = zeros(3, L);
    der_current_density_vec = zeros(3, L);
    E_mat = zeros(3, L, 2);
    E_mat(:,:,1) = E_1(:,:);
    E_mat(:,:,2) = E_2(:,:);
    for t = 1:L
        for i=1:2
            for j = 1:2
                for k = 1:2
                    p(:, t) = p(:, t) + eps0 .* chi_3 .* dot(E_mat(:,t,i), E_mat(:,t,j))*E_mat(:,t,k);
                end
            end
        end
    end
    for l=1:3
        current_dens_vec(l,:) = cent_diff_n(p(l,:), delta_t, 3);
        der_current_density_vec(l,:) = cent_diff_n(current_dens_vec(l,:) , delta_t, 3);
    end
    der_j_kerr = der_current_density_vec;
end