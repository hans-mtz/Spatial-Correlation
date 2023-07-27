function [beta_hat, u_hat] = ols(y,X,Z)
    %---------OLS-------------%
    T = length(y);
    % XX = X'*X/T;
    % Xy = X'*y/T;
    % XX_inv = inv(XX);
    % beta_hat = XX_inv*Xy;

    ZX = Z'*X/T;
    Zy = Z'*y/T;
    ZX_inv = inv(ZX);
    
%     beta_hat = ZX\Zy; %Matlab says this is faster than ZX_inv*Zy;
    beta_hat = ZX_inv*Zy;
    u_hat = y - X*beta_hat;
    % fprintf('Checking moment condition E[Zu]=0 %5.4g \n', Z'*u_hat/T);
    % fprintf('Checking moment condition E[Xu]=0 %5.4g \n', X'*u_hat/T);
end