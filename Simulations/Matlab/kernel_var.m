function [SE]= kernel_var(u_hat,X,Z,D_mat,L)
[T, k] = size(X);
ZX = Z'*X/T;
ZX_inv = inv(ZX);
% ZX_inv = ZX\1;
%------Variance kernel estimator------------%
Zu_hat = Z.*repmat(u_hat,1,k);
% Xu_hat = X.*repmat(u_hat,1,k);

V_hat = zeros(k); % Initializing V matrix

for i=1:T
    for j=find(D_mat(i,:)<=L)
%       for j=i
        V_hat = V_hat + Zu_hat(i,:)'*Zu_hat(j,:);
    end
end

V_hat0 = V_hat/T;

% V_kernel = (XX_inv)*V_hat*(XX_inv)'*T/(T-k); Stata's DF
V_kernel = (ZX_inv)*V_hat0*(ZX_inv)'/T;

SE = sqrt(diag(V_kernel));

end
