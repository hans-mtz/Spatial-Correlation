function bic_value = get_opt_splines(q_in,s,o,y,X,cholesky_flag)
    % Matrix of BSplines on locations
    q= round(q_in);
    [S, ~, ~]=get_bspline_mat(s,q,o,1);
    X_bs = [X(:,2:end) S];
    [~, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag);
    bic_value=bic(u_hat_bs,X_bs,"hansen");
end