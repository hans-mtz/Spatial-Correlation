function [SE]= kernel_var(u_hat,X,Z,D_mat,L,cholesky_flag)
[T, k] = size(X);
ZX = Z'*X/T;

switch cholesky_flag
    case 'ldl'
        
        [LA, DA, PA] = ldl(ZX);
        ZX_inv = PA*(LA'\(DA\(LA\(PA'*eye(k)))));
        
    case 'chol'
        R = chol(ZX);
        ZX_inv = R\(R'\eye(k));
    otherwise
        ZX_inv = inv(ZX);
end
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
V_kernel = (ZX_inv)*V_hat0*(ZX_inv)'/(T-k); % Adjusting for d.f. (before only T)

SE = sqrt(diag(V_kernel));

end
