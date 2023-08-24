function [n_splines, min_bic] = get_opt_splines_grid(s,o,y,X,cholesky_flag)
    bic_vec = [];
    grid = 5:5:90;
    for i=grid
        [S, ~, ~]=get_bspline_mat(s,i,o,1);
        X_bs = [X(:,2:end) S];
        [~, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag);
        bic_vec=[bic_vec bic(u_hat_bs,X_bs,"hansen")];
    end
    [min_bic, index]=min(bic_vec);
    n_splines = grid(index);

end