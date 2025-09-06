function sigma_tau = get_sigma_tau(tau,D_mat,type)
    arguments
        tau (2,1) double % tau parameters for the covariance matrix
        D_mat (:,:) double % Distance matrix
        type char {mustBeMember(type,{'dgp','exp','dsqr'})} = 'dsqr' % Type of covariance function
    end
    
    switch type
        case 'dgp'
            n = size(D_mat, 1);
            sigma_tau = exp(tau(1).*eye(n)).*exp(-tau(2).*D_mat);
        
        case 'exp'
            sigma_tau = exp(tau(1)).*exp(-tau(2).*D_mat);

        otherwise % 'dsqr'
            sigma_tau = exp(tau(1)).*exp(-tau(2).*(D_mat.^2));
    end
end