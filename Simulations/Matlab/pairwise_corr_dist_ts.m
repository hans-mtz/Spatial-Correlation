%% Defining groups based on distance %%
clear;

rng(333);

L= [0.01 0.05 0.1 0.5 0.8];
T=250;
beta = [1; 1.5];

%% Locations %%

s = (1:T)'./T;

D = getdistmat(s,0);

% generate data      
[y, X, u] = DGP_ts(beta, 0.79,T,0);
% Running OLS
[beta_hat, u_hat] = ols(y,X,X,'chol');

S = repmat(u_hat,1,T);

%% Pairwise correlation %%
rho_bar_dist = NaN(length(L),2);
k=1;
for l=L
rho=[];%NaN(T,T);
rho_no_over=[];%NaN(T,T);

for i=1:T
    for j=1:T
        A=S(D(i,:)<=l,i);
        B=S(D(j,:)<=l,j);
        
        if length(A)>length(B)
            B_temp = B;
            B = NaN(length(A),1);
            B(1:length(B_temp)) = B_temp;
        elseif length(A)<length(B)
            A_temp = A;
            A = NaN(length(B),1);
            A(1:length(A_temp)) = A_temp;            
        end
        
        if (isempty(intersect(A,B))) && (length(A)==length(B))     
            rho_no_over = [rho_no_over corr(A,B,'Rows','all')];
        end
        
        if (j~=i) && (length(A)==length(B))     
            rho = [rho corr(A,B,'Rows','all')];
        end
    end
end

rho_bar_dist(k,1)=mean(rho,'all','omitnan');
rho_bar_dist(k,2)=mean(rho_no_over,'all','omitnan');
k=k+1;
end
rho_bar_dist

tbl = cat(2,L',rho_bar_dist(:,1));

% 
writetable(array2table(tbl),'../Products/pairwise_corr.csv');

