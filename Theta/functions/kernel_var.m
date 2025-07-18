function [SE]= kernel_var(u_hat,X,Z,D_mat,L,s,kernel,cholesky_flag, fix)


    if L <= 0
        SE = HR_var(u_hat,X,Z);
        fprintf('Using HR variance estimator as L <= 0\n');
        return;
    end

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

    switch kernel
        case 'uniform'
            for i=1:T
                for j=find(D_mat(i,:)<=L)
            %       for j=i
                    V_hat = V_hat + Zu_hat(i,:)'*Zu_hat(j,:);
                end
            end
        case 'gaussian'
            K = normpdf(D_mat,0,L/2)./normpdf(0,0,L/2);
            for i=1:T
                for j=1:T
                    V_hat = V_hat + Zu_hat(i,:)'*Zu_hat(j,:)*K(i,j);
                end
            end
        otherwise %'triangle'
            K = triangle_weights(s,L);
            for i=1:T
                for j=find(K(i,:)> 0)
                    V_hat = V_hat + Zu_hat(i,:)'*Zu_hat(j,:).*K(i,j);
                end
            end
    end

    V_hat0 = V_hat/T;

    % V_kernel = (XX_inv)*V_hat*(XX_inv)'*T/(T-k); Stata's DF
    V_kernel = (ZX_inv)*V_hat0*(ZX_inv)'/(T-k); % Adjusting for d.f. (before only T)

    if (fix)
        %vcov fix if negative elements in the diagonal a la Cameron, Gelbech and Miller (2011)
        V_k = V_plus(V_kernel);
    else
        V_k = V_kernel;
    end

    SE = sqrt(diag(V_k));

end

function V_hat_plus= V_plus(V_hat)

    if any(diag(V_hat)<0)
        [V, D] = eig(V_hat);
        D_plus = diag(max(D));
        V_hat_plus = V*D_plus*V';
    else
        V_hat_plus = V_hat;
    end

end
