function [index, min_value] = get_opt_n_pcs_grid(y,X,W,D_mat,method,cholesky_flag)
    value_vec = NaN(size(W,2),1);
    % betas = NaN(size(W,2),1);
    % residuals = NaN(size(W,2),1);

    % grid = 1:size(W,2);

    switch method
        case 'bic'

            for i=1:size(W,2)
                % [W]=get_PC(i,s,q_max,o);
                S = W(:,1:i);
                X_bs = [X S];
                [beta_hat_bs, u_hat_bs] = ols(y,X_bs,X_bs,cholesky_flag);
                value_vec(i) = bic(u_hat_bs,X_bs,"hansen");
                % betas(i) = beta_hat_bs;
                % residuals(i) = u_hat_bs;
            end
            
        otherwise

            for i=1:size(W,2)
                S = W(:,1:i);
                X_bs = [X S];
                [beta_hat_bs, u_hat_bs] = ols(y,X_bs,X_bs,cholesky_flag);

                nn_cor = get_nn_corr(u_hat_bs,D_mat,1);
                value_vec(i) = abs(nn_cor-0.0); % 0.05 is the target correlation
                % betas(i) = beta_hat_bs;
                % residuals(i) = u_hat_bs;
            end

    end

    [min_value, index] = min(value_vec);
    % q = index;
    % beta_hat = betas(index);
    % u_hat = residuals(index);

end