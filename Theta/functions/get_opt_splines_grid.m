function [n_splines, min_value] = get_opt_splines_grid(s,o,y,X,method,cholesky_flag)
    % Returns the optimum number of splines that minimizes the criterion
    % provided by 'method' given a vector of locations s, the order of splines o,
    % the output vector y, the input matrix X. The cholesky_flag chooses how the OLS
    % function inverts the X'X matrix

    n_s=size(s,2); %number of locations per observation
    D_mat = getdistmat(s,false); % Distance matrix
    if n_s>1
        grid=4:2:12;
    else
        grid = 5:5:90;
    end
    value_vec = NaN(size(grid));
    j=1;
    switch method
        case 'bic'

            for i=grid
                [S, ~, ~]=get_bsplines(s,i,o);
                X_bs = [X S(:,1:end-1)];
                [~, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag);
                value_vec(j)=bic(u_hat_bs,X_bs,"hansen");
                j=j+1;
            end
            
        otherwise

            for i=grid
                [S, ~, ~]=get_bsplines(s,i,o);
                X_bs = [X S(:,1:end-1)];
                [~, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag);
                if (n_s>1)
                    gamma=get_nn_corr(u_hat_bs,D_mat,1);
                else
                    [u_t, u_t_1] = sort_u(u_hat_bs,s);
                    beta_ar1 = ols(u_t, u_t_1, u_t_1, cholesky_flag);
                    gamma=beta_ar1(2);
                end
                value_vec(j)=abs(gamma-0.05);
                j=j+1;
            end

    end

    [min_value, index]=min(value_vec);
    n_splines = grid(index);

end