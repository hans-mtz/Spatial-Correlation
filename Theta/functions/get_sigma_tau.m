function sigma_tau = get_sigma_tau(tau,D_mat)
    n = size(D_mat, 1);
    % sigma_tau = exp(tau(1)).*exp(-tau(2).*D_mat);
    sigma_tau = exp(tau(1)).*exp(-tau(2).*(D_mat.^2));
    % sigma_tau = exp(tau(1).*eye(n)).*exp(-tau(2).*D_mat);
end