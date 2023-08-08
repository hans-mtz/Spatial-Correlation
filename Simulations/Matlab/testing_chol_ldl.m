% Setting up parameters
global T k beta rho_bar L
T = 250;
k= 3;
beta = [ 1; 1.5; 0.5];
rho_bar = 0.03; % Average pairwise correlation

% Fixing locations
s = rand(T,1); % vector of locations
S=get_bspline_mat(s,4,2); % Matrix of BSplines on locations

% generate data      
[y, X, D_mat] = DGP(beta,k,s,rho_bar,2);

%cholesky factorization and inversion

A = X'*X/T;
Xy = X'*y/T;
R = chol(A);
L_chol = R';
b = R'\Xy;
beta_hat_chol = R\b;
beta_hat_chol = R\(R'\Xy);
beta_hat_chol = ols(y,X,X, true);
beta_hat_ols = ols(y,X,X, false);

[LA, DA] = ldl(A);
beta_hat_ldl = LA'\(DA\(LA\Xy));
[LA, DA, PA] = ldl(A);
beta_hat_ldl_p = PA*(LA'\(DA\(LA\(PA'*Xy))));
%% Cholesky efficient %%

A = X'*X/T;
Xy = X'*y/T;
R = chol(A);
L_chol = R'; %chol(A,'lower');
S_temp = zeros(size(A));
for i=1:length(L_chol)
    for j=i
        S_temp(i,j)= 1/L_chol(i,i);
    end
end

A_inv = NaN(size(A));

for i=1:length(A_inv)
    A_inv(:,i) = R\S_temp(:,i);
end

beta_hat_chol_eff = A_inv*Xy;


        

