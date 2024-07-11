function SE_HR = HR_var(u_hat, X, Z)   
% HR_var - Estimating Heteroscedastic Robust Variance
% using Stata's DF adjustment
%
% Syntax: SE = HR_var(u_hat, X, )   
%
% Long description
[T, k] = size(X);
ZX = Z'*X/T;
ZX_inv = inv(ZX);
% ZX_inv = ZX\1;
%------HR Variance estimator------------%
Zu_hat = Z.*repmat(u_hat,1,k)/T;
V_hat = Zu_hat'*Zu_hat;
V_HR = ZX_inv*V_hat*ZX_inv';
V_HR = V_HR*T/(T-k);
SE_HR = sqrt(diag(V_HR));

    
end