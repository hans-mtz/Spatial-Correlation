function [statistic, se] = get_statistic(beta_hat,beta0,u_hat,X,D_mat,cutoff,s)
    se = kernel_var(u_hat,X,X,D_mat,cutoff,s,'gaussian','chol',false);

    H_0 = beta_hat - beta0;
    if any(se~=real(se))
        warning('Some standard errors are imaginary. Returning NaN for statistic.');
        statistic = NaN(size(H_0));
    else
        statistic = H_0 ./ se;
    end
    statistic = statistic(:); % Ensure it is a column vector
end