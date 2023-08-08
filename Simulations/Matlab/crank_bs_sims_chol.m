% Simulations - Rejection Frequencies for Kernel Var Estimator%
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


% Setting up parameters
global T k beta rho_bar L
T = 250; % obs
beta = [ 1; 1.5];
k= length(beta); % number of covariates
rho_bar = 0.03; % Average pairwise correlation
L = [0.015 0.03 0.05]; % Vector of cutoff distances
R = [4 12 24 36];% 90 120 180]; % Number of triangles
% R = 40:48; % 43 for order 1 splines
% R = 64:84; % 82 for order 2 splines
cholesky_flag = 'ldl';
N_order_splines = 2;
N_models = 2;
N_reps = 250;

N_cutoffs = length(L);
N_triangles = length(R);
N_estimators = 2+N_cutoffs+N_cutoffs*N_triangles*N_order_splines; 


% Fixing locations
s = rand(T,1); % vector of locations
rej_freq = NaN(k,N_estimators,N_models);

tic
for m=2:3
    
    rej_freq_m_i = zeros(N_reps,k,N_estimators);

    for r=1:N_reps
        fprintf('model %d, rep %d \n',m,r);
        % generate data      
        [y, X, D_mat] = DGP(beta,s,rho_bar,m);
%         writetable(table(y, X, s), strcat('../Stata/data_',num2str(m),'.csv')); % Saving
        
        % Running OLS
        [beta_hat, u_hat] = ols(y,X,X,cholesky_flag);

        % Getting SE from HR and SCPC estimatots
        se_hr = HR_var(u_hat,X,X);
        [se_scpc, cv_scpc] = scpc_var(u_hat,beta_hat,s,X,rho_bar,0.95,0);
        
        % Testing HR and SCPC
        H_0 = (beta_hat - beta);
        t_stat_hr = H_0./se_hr;
        rej_hr = abs(t_stat_hr) > 1.96;
        rej_freq_m_i(r,:,1) = rej_hr;

        t_stat_scpc = H_0./se_scpc;
        rej_scpc = abs(t_stat_scpc) > abs(cv_scpc);
        rej_freq_m_i(r,:,2) = rej_scpc; 

        % Kernel for l\in L
        i=3;
        for l=L
            se_kernel = kernel_var(u_hat,X,X,D_mat,l, cholesky_flag);
            t_stat_kernel = H_0./se_kernel;
            rej_kernel = abs(t_stat_kernel) > 1.96;
            rej_freq_m_i(r,:,i) = rej_kernel;
            i=i+1;
        end
        % Kernel + BSplines
        j=i; %2+length(L)+1;
        for o=1:N_order_splines
            for q=R %Number of splines
                for p=L %Cutoff lengths for Kernel
%                     disp(q);
                    S=get_bspline_mat(s,q,o); % Matrix of BSplines on locations
                    X_bs = [X(:,2:end) S];
                    [beta_hat_bs, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag); % OLS + BSplines, NO Constant
                    se_kernel_bs = kernel_var(u_hat_bs,X_bs,X_bs,D_mat,p, cholesky_flag);
                    t_stat_kernel_bs = (beta_hat_bs(1:k-1)-beta(2:k))./se_kernel_bs(1:k-1);
                    rej_kernel_bs = abs(t_stat_kernel_bs) > 1.96;
                    rej_freq_m_i(r,:,j) = [NaN rej_kernel_bs];
                    j=j+1;
                    
                end
            end
        end 
    end
    rej_freq(:,:,m-1) = sum(rej_freq_m_i,1)./N_reps;
end
toc

    
results_table = array2table(...
    cat(2,rej_freq(:,:,1)',rej_freq(:,:,2)'));
results_table

% results_table = array2table(rej_freq(:,:,1)');
% results_table
writetable(results_table,'../Products/crank_bs_sims_chol.csv');
% writetable(table(y, X, s), '../Products/example_matlab.csv'); % Saving
        