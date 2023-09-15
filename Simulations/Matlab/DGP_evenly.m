function [y,X,u] = DGP_evenly(beta,coef_vec,D_mat,L)
    T=size(D_mat,1);
    burn = 100;
    X_new = zeros(T,1);
    u_new = zeros(T,1);
    epsilon = randn(T,2,5+burn);
    X_old = epsilon(:,1,1);
    u_old = epsilon(:,2,1);
    
    if L>0 && L<=1
        for r=2:(4+burn)
            
            for i=1:T
                coef_j = coef_vec;
%                 neighbours = D_mat(i,D_mat(i,:)<=L);
%                 sorted_neigs = sort(neighbours);
%                 vals = D_mat(i,1:end ~= i)<=L;
                
                for j=find(D_mat(i,:)<=L)
                    if j~=i
                        X_new(i) = X_new(i)+X_old(j)*coef_j;
                        u_new(i) = u_new(i)+u_old(j)*coef_j;
                        coef_j=coef_j*coef_vec;
                    end
                end
            end
            
            X_new=X_new+epsilon(:,1,r-1);
            u_new=u_new+epsilon(:,2,r-1);
            X_old=epsilon(:,1,r);
            u_old=epsilon(:,2,r);
        end
    else 
        for r=2:(4+burn)
            
            for i=1:T
                
                val = min(D_mat(i,1:end ~= i));
                coef_j = coef_vec;
                
                for j=find(D_mat(i,:)<=val)
                    if j~=i
                        X_new(i) = X_new(i)+X_old(j)*coef_j;
                        u_new(i) = u_new(i)+u_old(j)*coef_j;
                        % coef_j=coef_j*coef_vec;
                    end
                end
            end
            
            X_old=X_new;
            u_old=u_new;
            X_new=epsilon(:,1,r);
            u_new=epsilon(:,2,r);
        end
        
    end
    X=[ones(T,1) X_old];
    u=u_old;
    y=X*beta+u;
end