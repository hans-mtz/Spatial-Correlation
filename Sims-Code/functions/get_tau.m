function [tau, beta_hat, exitflag] = get_tau(y,X,D_mat,verbose)
    arguments
        y (:,1) double
        X (:,:) double
        D_mat (:,:) double
        verbose logical = false
    end

    % Get residuals
    [beta_hat, u_hat] = ols(y, X, X,'chol');
    % Get tau
    fun = @(tau) llh(tau, u_hat, D_mat);
    options = optimset('Display', 'notify', 'TolFun', 1e-4, 'TolX', 1e-4,'MaxFunEvals', 1000, 'MaxIter', 1000); %'Display', 'off'
    tau0 = [0.1; 1]; % Initial guess for tau
    [tau1, ~, ~] = fminsearch(fun, tau0, options);
    if verbose
        fprintf('Initial tau: %f, %f\n', tau1(1), tau1(2)); 
    end
    options = optimset('Display', 'notify'); %'Display', 'notify' 'iter'
    [tau, fval, exitflag] = fminsearch(fun, tau1, options);
    if verbose
        % fprintf('Optimized tau: %f, %f\n', tau(1), tau(2));
        if exitflag <= 0 
            fprintf('tau: %f, %f\n', tau(1), tau(2));
            warning('Optimization did not converge: %s', optimget(options, 'Display'));
        else
            fprintf('Optimization converged with function value: %f\n', fval);
            fprintf('tau: %f, %f\n', tau(1), tau(2));
        end
    end
    % tau = [exp(tau(1)), tau(2)]; % Return tau in the desired format
end

