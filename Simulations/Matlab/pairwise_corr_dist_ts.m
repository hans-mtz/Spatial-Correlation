%% Defining groups based on distance %%
clear;

rng(333);

L= [0.015 0.05 0.1 0.5 0.9 1];
T=250;
beta = [1; 1.5];

%% Locations %%

s = (1:T)'./T;

D = getdistmat(s,0);

% generate data      
[y, X, u] = DGP_ts(beta, 0.79,T,0);
% Running OLS
[beta_hat, u_hat] = ols(y,X,X,'chol');
% u_hat = u;

% S = repmat(u_hat,1,T);

%% Pairwise correlation %%

rho_bar = NaN(length(L),1);
k=1;
for l=L
rho_i=NaN(T,1);
for i=1:T
    rho=[];%NaN(T,T);
    index=find(D(i,:)<=l);
    mu = mean(u_hat(index));
    sigma = var(u_hat(index));

    for j=index
        if i~=j
            rho = [rho (u_hat(i)*u_hat(j)-mu)];
            
        end
    end
%     rho_i(i) = (mean(rho))/(sigma);
%     rho_i(i) = (mean(rho)-(mu*mu))/norm(u_hat(index));
    rho_i(i) = mean(rho);
end
rho_bar(k) = sum(rho_i)./((T-1));
k=k+1;
end
rho_bar

tbl = cat(2,L',rho_bar)

%% Saving %%

writetable(array2table(tbl),'../Products/pairwise_corr.csv');

