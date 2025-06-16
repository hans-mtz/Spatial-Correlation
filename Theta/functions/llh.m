function [val] = llh(tau, u_hat, D_mat)
    % arguments
    %     tau (1,2) double % Parameters to optimize
    %     u_hat (:,1) double % Residuals from OLS regression
    %     D_mat (:,:) double % Distance matrix
    % end
    n = length(u_hat);
    % sigma = (tau(1)-1).*eye(n) + exp(-tau(2).*D_mat);
    sigma = exp(tau(1).*eye(n)).*exp(-tau(2).*D_mat);
    [R, flag ]= chol(sigma); % Cholesky decomposition
    % Ensure sigma is positive definite
    if (flag ~= 0)
        val=Inf;
        return; % Return negative infinity if not positive definite
    end
    val = -0.5.*log(det(sigma)) - 0.5.*u_hat'*(R\(R'\u_hat));
    val = -val; % Negate for minimization
    % Note: The function is negated because fminsearch minimizes the function.
end