function [tau, beta_hat] = get_tau(y,X,D_mat)
    arguments
        y (:,1) double
        X (:,:) double
        D_mat (:,:) double
    end

    % Get residuals
    [beta_hat, u_hat] = ols(y, X, X,'chol');
    % Get tau
    fun = @(tau) llh(tau, u_hat, D_mat);
    options = optimset('Display', 'off', 'TolFun', 1e-8, 'TolX', 1e-8,'MaxFunEvals', 1000, 'MaxIter', 1000);
    tau0 = [1; 1]; % Initial guess for tau
    [tau1, ~, ~] = fminsearch(fun, tau0, options);
    options = optimset('Display', 'notify');
    [tau, fval, exitflag] = fminsearch(fun, tau1, options);
    if exitflag <= 0
        warning('Optimization did not converge: %s', optimget(options, 'Display'));
    else
        fprintf('Optimization converged with function value: %f\n', fval);
    end
    % tau = [exp(tau(1)), tau(2)]; % Return tau in the desired format
end

