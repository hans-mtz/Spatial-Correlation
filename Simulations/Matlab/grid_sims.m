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
rng(333) % Setting seed
excercise = 'dd_2x';

%% Setting up parameters
global T k beta rho_bar L
N_reps = 1000;
T = 250; % obs
beta = [ 1; 1.5];
k= length(beta); % number of covariates
rho_bar = 0.03; % Average pairwise correlation
L = [0.015 0.03 0.05]; % Vector of cutoff distances
% R = [4 12 24 36];% 90 120 180]; % Number of triangles
R = [4 8 12:12:72];%[4 12 24 36 48 64 88 102 124];
% R = 40:48; % 43 for order 1 splines
% R = 64:84; %82 for order 2 splines
% models = ['iid' 'AR(1)'];
cholesky_flag = 'chol';
spatial = 0; % 1 s ~ U(0,1); 0 s ~ [1:T]./T
N_order_splines = 2;
N_cutoffs = length(L);
N_triangles = length(R);
N_estimators = (1+N_cutoffs)+(2+N_cutoffs)*N_triangles*N_order_splines; 
N_ols = 1+N_triangles*N_order_splines;

%% Setting locations -------
if spatial==1
    % Fixing locations
    s = rand(T,1); % vector of locations
else % Time series
    s = (1:T)'./T;
end

rej_freq = NaN(k,N_estimators);
avg_e_ar1 = NaN(4,N_ols);
rej_freq_m_i = NaN(N_reps,k,N_estimators);
e_ar1_betas = NaN(N_reps,4,N_ols);

tic
for r=1:N_reps
    if N_reps<100 || mod(r,100)==0
        fprintf('rep %d \n',r);
    end
    % generate data      
%         [y, X, D_mat] = DGP(beta,s,rho_bar,m,spatial); %y and X change every draw
    [y, X] = DGP_ts(beta, 0.79,T,0);
%      writetable(table(y, X, s), strcat('../Products/data_ts.csv')); % Saving
    
    % Running OLS
    [beta_hat, u_hat] = ols(y,X,X,cholesky_flag);
    
    % Testing HR and SCPC
    H_0 = (beta_hat - beta);
    
    % Residuals' AR(1)
    [u_t, u_t_1] = sort_u(u_hat,s);
    e_ar1_betas(r,1:2,1) = ols(u_t, u_t_1, u_t_1, cholesky_flag);
    e_ar1_betas(r,3,1) = bic(u_hat,X,"hansen");
    e_ar1_betas(r,4,1) = bic(u_hat,X,"damian");
    % Getting SE from HR and SCPC estimatots
    se_hr = HR_var(u_hat,X,X);
    t_stat_hr = H_0./se_hr;
    rej_hr = abs(t_stat_hr) > 1.96;
    rej_freq_m_i(r,:,1) = rej_hr;
    
%         [se_scpc, cv_scpc] = scpc_var(u_hat,beta_hat,s,X,rho_bar,0.95,0);
%         t_stat_scpc = H_0./se_scpc;
%         rej_scpc = abs(t_stat_scpc) > abs(cv_scpc);
%         rej_freq_m_i(r,:,2) = rej_scpc; 

    % Kernel for l\in L
    D_mat = getdistmat(s,false);
    i=2;
    for l=L
        se_kernel = kernel_var(u_hat,X,X,D_mat,l,cholesky_flag);
        t_stat_kernel = H_0./se_kernel;
        rej_kernel = abs(t_stat_kernel) > 1.96;
        rej_freq_m_i(r,:,i) = rej_kernel;
        i=i+1;
    end
    % Kernel + BSplines
    j=i; %2+length(L)+1;
    n=2;
    if r==1
        cutoff_dist_vec = [];
    end
    for o=1:N_order_splines
        for q=R %Number of splines

            [S, delta]=get_bspline_mat(s,q,o,0); % Matrix of BSplines on locations
            if r==1
                cutoff_dist_vec=[cutoff_dist_vec; NaN; L'; 2*delta*o];
            end
            X_bs = [X(:,2:end) S]; %Taking constant off, adding splines to regressors
            [beta_hat_bs, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag); % OLS + BSplines, NO Constant
            %Null hypothesis against the truth
            H_0_bs = beta_hat_bs(1:k-1)-beta(2:k);
            % Getting SE from HR and SCPC estimatots
            se_hr_bs = HR_var(u_hat_bs,X_bs,X_bs);
            t_stat_hr_bs = H_0_bs./se_hr_bs(1:k-1);
            rej_hr_bs = abs(t_stat_hr_bs) > 1.96;
            rej_freq_m_i(r,:,j) = [NaN rej_hr_bs];
            % Residuals' AR(1)
            [u_t_bs, u_t_1_bs] = sort_u(u_hat_bs,s);
            e_ar1_betas(r,1:2,n) = ols(u_t_bs, u_t_1_bs, u_t_1_bs, cholesky_flag);
            e_ar1_betas(r,3,n) = bic(u_hat_bs,X_bs,"hansen");
            e_ar1_betas(r,4,n) = bic(u_hat_bs,X_bs,"damian");
            j=j+1;
            for p=[L 2*delta*o]%Cutoff lengths for Kernel
                
                se_kernel_bs = kernel_var(u_hat_bs,X_bs,X_bs,D_mat,p, cholesky_flag);
                t_stat_kernel_bs = (H_0_bs)./se_kernel_bs(1:k-1);
                rej_kernel_bs = abs(t_stat_kernel_bs) > 1.96;
                rej_freq_m_i(r,:,j) = [NaN; rej_kernel_bs];
                j=j+1;
                
            end
            % advancing error AR(1) beta matrix index
            n=n+1;
        end
    end 
end

rej_freq(:,:) = sum(rej_freq_m_i,1)./N_reps;
avg_e_ar1(:,:) = sum(e_ar1_betas,1)./N_reps;
cutoff_dist_vec
toc

%% Preparing tables %%
    
results_table = array2table(...
    rej_freq(:,:)',...
    'VariableNames', {'Cons.' 'Beta'});
results_table

ear1_table = array2table(...
    avg_e_ar1(:,:)',...
    'VariableNames', {'Cons.' 'Gamma' 'BIC_h' 'BIC_d'});
ear1_table

%% Saving results %%
save(['grid_' excercise '.mat'])

writetable(results_table,['../Products/grid_sims_res_' excercise '.csv']);
writetable(ear1_table,['../Products/grid_sims_ar1_' excercise '.csv']);
writetable(array2table(cutoff_dist_vec) ,['../Products/grid_sims_dist_v_' excercise '.csv']);

