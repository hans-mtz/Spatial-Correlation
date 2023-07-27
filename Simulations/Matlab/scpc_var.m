function [se, cv] = scpc_var(u_hat, beta_hat, s, X,rhobar,ci_level,latlongflag)
    [n, k] = size(X);
    Xu = X.*repmat(u_hat,1,k);
    XX = (X'*X);
    S_XX = XX/n;
    S_XX_inv = inv(S_XX);
    y_data = repmat(beta_hat',n,1)+Xu*S_XX_inv';

    rslt = scpc(y_data,s,latlongflag,rhobar,ci_level);  % Run SCPC
    se = rslt.se_beta_hat;
    cv = rslt.cv;

end