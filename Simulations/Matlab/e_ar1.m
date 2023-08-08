%% Simulations - Rejection Frequencies for Kernel Var Estimator%
clear;

if batchStartupOptionUsed
    addpath(genpath('./SCPC_Matlab_Implementation_20211123'))
end
%--------------Models-------------------%
% Model 1 
% e~G_exp(c_rho_bar) X~[1 N(0,1) N(0,1)]
% Model 2 
% e~G_exp(c_rho_bar) 
% X~[1 G_exp(c_rho_bar) G_exp(c_rho_bar)]
% Model 3
% e~G_exp(c_rho_bar) 
% X~[1 G_exp(c_rho_bar) G_exp(c_rho_bar)]
% Y demeaned in estimation

%--------------Estimators---------------%
%HR, Kernel, Kernel + splines, SCPC
rng(734) % Setting seed


%% Setting up parameters
global T k beta rho_bar
T = 250; % obs
beta = [ 1; 1.5; 0.5];
k= length(beta); % number of covariates
rho_bar = 0.03; % Average pairwise correlation
R = [4 12 24 36];% 90 120 180]; % Number of triangles
% R = [64 32 12 5];
cholesky_flag = 'chol';
N_order_splines = 2;

N_splines = length(R);
N_estimators = 1+N_splines*N_order_splines; 


%% Fixing locations
s = rand(T,1); % vector of locations

beta_hats = NaN(k,N_estimators);
u_ar1_beta = NaN(2,N_estimators);

tic

[y, X, D_mat] = DGP(beta,k,s,rho_bar,2);

[beta_hat, u_hat] = ols(y,X,X,cholesky_flag);

beta_hats(:,1) = beta_hat;

[u_t, u_t_1] = sort_u(u_hat,s);

[beta_u_hat] = ols(u_t, u_t_1, u_t_1, cholesky_flag);

u_ar1_beta(:,1) = beta_u_hat;
%% Running loop

j=2; %2+length(L)+1;
for o=1:N_order_splines
    for q=R %Number of splines

            
        S=get_bspline_mat(s,q,o); % Matrix of BSplines on locations

        [beta_hat_bs, u_hat_bs] = ols(y,[X(:,2:end) S],[X(:,2:end) S], cholesky_flag); % OLS + BSplines, NO Constant

        beta_hats(2:end,j) = beta_hat_bs(1:k-1);

        [u_t, u_t_1] = sort_u(u_hat_bs,s);

        [beta_u_hat_bs] = ols(u_t, u_t_1, u_t_1, cholesky_flag);

        u_ar1_beta(:,j) = beta_u_hat_bs;

        j=j+1;
                

    end
end 

toc

%% Saving results
    
results_table = array2table(...
    cat(2,beta_hats',u_ar1_beta'),...
    'VariableNames',{'Cons' 'Beta1' 'Beta2' 'AR1 Cons' 'AR1 Beta'});
results_table

writetable(results_table,'../Products/e_ar1.csv');
