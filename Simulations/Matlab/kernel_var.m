function [beta_hat, SE] = kernel_var(y,X,D_mat,k,L,T)
%---------OLS-------------%
XX = X'*X;
% ZX = Z'*X;
% Zy = Z'*y;
Xy = X'*y;
% ZX_inv = inv(ZX);
XX_inv = inv(XX);

% beta_IV_hat = ZX_inv*Zy;
beta_hat = XX_inv*Xy;
u_hat = y - X*beta_hat;
% fprintf('Checking moment condition E[Zu]=0 %5.4g \n', Z'*u_hat/T);
% fprintf('Checking moment condition E[Xu]=0 %5.4g \n', X'*u_hat/T);


%------Variance kernel estimator------------%
% Zu_hat = Z.*repmat(u_hat,1,k);
Xu_hat = X.*repmat(u_hat,1,k);

V_hat = zeros(k); % Initializing V matrix

for i=1:T
    for j=find(D_mat(i,:)<=L)
%       for j=i
        V_hat = V_hat + Xu_hat(i,:)'*Xu_hat(j,:);
    end
end

% V_hat0 = V_hat/T; No need to make comparable to Stata's HR SE 
% when K=1 if i=j; 0 otherwise

% V_kernel = (XX_inv)*V_hat*(XX_inv)'*T/(T-k); Stata's DF
V_kernel = (XX_inv)*V_hat*(XX_inv)'/T;
SE = sqrt(diag(V_kernel));

% rt.beta_hat = beta_hat;
% rt.SE = SE;

end
