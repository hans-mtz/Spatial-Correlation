function [n_splines, min_value] = get_opt_splines_grid(s,o,y,X,method,cholesky_flag)
    value_vec = [];
    n_s=size(s,2);
    D_mat = getdistmat(s,false);
    if n_s>1
        grid=4:2:12;
    else
        grid = 5:5:90;
    end
    switch method
        case 'bic'

            for i=grid
                [S, ~, ~]=get_bsplines(s,i,o);
                X_bs = [X(:,2:end) S];
                [~, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag);
                value_vec=[value_vec bic(u_hat_bs,X_bs,"hansen")];
            end
            
        otherwise

            for i=grid
                [S, ~, ~]=get_bsplines(s,i,o);
                X_bs = [X(:,2:end) S];
                [~, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag);
                if n_s>1
                    gamma=get_nn_corr(u_hat_bs,D_mat,1);
                else
                    [u_t, u_t_1] = sort_u(u_hat_bs,s);
                    beta_ar1 = ols(u_t, u_t_1, u_t_1, cholesky_flag);
                    gamma=beta_ar1(2);
                end
                value_vec=[value_vec abs(gamma-0.05)];
            end

    end

    [min_value, index]=min(value_vec);
    n_splines = grid(index);

end