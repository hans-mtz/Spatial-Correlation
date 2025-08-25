function [beta_hat, u_hat] = ols(y,X,Z,chol_flag)
    %---------OLS-------------%
    T = length(y);
    % X = [ones(T,1) X_in]; % Adding intercept
    % Z = [ones(T,1) Z_in]; % Adding intercept
    % XX = X'*X/T;
    % Xy = X'*y/T;
    % XX_inv = inv(XX);
    % beta_hat = XX_inv*Xy;

    ZX = Z'*X/T;
    Zy = Z'*y/T;
    
    switch chol_flag
        case 'ldl'
            
            [LA, DA, PA] = ldl(ZX);
            beta_hat = PA*(LA'\(DA\(LA\(PA'*Zy))));
            
        case 'chol'
            
            R = chol(ZX);
            beta_hat = R\(R'\Zy);
            
        otherwise
            
            beta_hat = ZX\Zy; %Matlab says this is faster than ZX_inv*Zy;
            % ZX_inv = inv(ZX);
            % beta_hat = ZX_inv*Zy;
    end
    
    u_hat = y - X*beta_hat;
    % fprintf('Checking moment condition E[Zu]=0 %5.4g \n', Z'*u_hat/T);
    % fprintf('Checking moment condition E[Xu]=0 %5.4g \n', X'*u_hat/T);
end