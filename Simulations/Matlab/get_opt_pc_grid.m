function [q, min_value] = get_opt_pc_grid(q_max,s,o,y,X,method,cholesky_flag)
    value_vec = [];
    grid = 4:4:q_max;
    switch method
        case 'bic'

            for i=grid
                [W]=get_PC(i,s,q_max,o,1);
                X_bs = [X(:,2:end) W];
                [~, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag);
                value_vec=[value_vec bic(u_hat_bs,X_bs,"hansen")];
            end
            
        otherwise

            for i=grid
                [W]=get_PC(i,s,q_max,o,1);
                X_bs = [X(:,2:end) W];
                [~, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag);
                [u_t, u_t_1] = sort_u(u_hat_bs,s);
                beta_ar1 = ols(u_t, u_t_1, u_t_1, cholesky_flag);
                value_vec=[value_vec abs(beta_ar1(2)-0.05)];
            end

    end

    [min_value, index]=min(value_vec);
    q = grid(index);

end