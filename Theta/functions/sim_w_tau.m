function [y_sim, eps_sim] = sim_w_tau(tau,D_mat,X,beta_cand,B)
    arguments
        tau (2,1) double % tau parameters for the covariance matrix
        D_mat (:,:) double % Distance matrix
        X (:,:) double % Design matrix
        beta_cand (1,:) double % Candidate beta coefficients
        B (1,1) double = 1 % Number of simulations
    end
    n = size(D_mat, 1);
    sigma = exp(tau(1).*eye(n)).*exp(-tau(2).*D_mat);

    eps_sim = mvnrnd(zeros(n, 1), sigma,B)'; % Generate multivariate normal random variables
    y_sim = X*beta_cand' + eps_sim; % Simulated dependent variable
end