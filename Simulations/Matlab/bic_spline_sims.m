%% Simulations - Rejection Frequencies for Kernel Var Estimator%
clear;

if batchStartupOptionUsed
    addpath(genpath('./SCPC_Matlab_Implementation_20211123'))
end
%--------------Models-------------------%
% Model 1 iid 
% e~G_exp(c_rho_bar) X~[1 N(0,1) N(0,1)]
% Model 2 Spatial
% e~G_exp(c_rho_bar) 
% X~[1 G_exp(c_rho_bar) G_exp(c_rho_bar)]
% Model 3 Spatial - Demeaned outcome
% e~G_exp(c_rho_bar) 
% X~[1 G_exp(c_rho_bar) G_exp(c_rho_bar)]
% Y demeaned in estimation

%--------------Estimators---------------%
%HR, Kernel, Kernel + splines, SCPC
rng(333) % Setting seed
excercise = 'dd_2x';%'grid';
grid = 0;

%% Setting up parameters
global T k beta rho_bar L
N_reps = 1000;
T = 250; % obs
beta = [ 1; 1.5];
k= length(beta); % number of covariates
rho_bar = 0.03; % Average pairwise correlation
L = [0.015 0.03 0.05]; % Vector of cutoff distances
% R = [4 12 24 36 48 64 88 102 124];% 90 120 180]; % Number of triangles
% R = 40:48; % 43 for order 1 splines
% R = 64:84; % 82 for order 2 splines
cholesky_flag = 'chol';
N_order_splines = 2;


N_cutoffs = length(L);
% N_triangles = length(R);
N_estimators = (3+N_cutoffs)+(2+N_cutoffs)*N_order_splines; 
N_ols = 1+N_order_splines;

% Fixing locations
s = rand(T,1); % vector of locations

%% Loop
rej_freq = NaN(k,N_estimators);
avg_e_ar1 = NaN(5,N_ols);
t_stat_array = NaN(N_reps,k, N_estimators);
num_array = NaN(N_reps,k, N_estimators);
den_array = NaN(N_reps,k, N_estimators);
rej_freq_m_i = NaN(N_reps,k,N_estimators);
e_ar1_betas = NaN(N_reps,5,N_ols);
splines_array = NaN(N_reps,N_order_splines+2);

tic
fprintf('Running simulations\n');
for r=1:N_reps
    
    % generate data      
    [y, X, D_mat] = DGP(beta,s,rho_bar,1);

%         writetable(table(y, X, s), strcat('../Stata/data_splines.csv')); % Saving
    for m=1%1:2

        if N_reps<100 || mod(r,100)==0
            fprintf('model %d, rep %d \n',m,r);
        end
        % Running OLS
        if m==2
            [beta_hat, u_hat] = ols(y-mean(y),X,X,cholesky_flag);
        else
            [beta_hat, u_hat] = ols(y,X,X,cholesky_flag);
        end

        % Residuals' AR(1)
        [u_t, u_t_1] = sort_u(u_hat,s);
        e_ar1_betas(r,1:2,1) = ols(u_t, u_t_1, u_t_1, cholesky_flag);
        e_ar1_betas(r,3,1) = bic(u_hat,X,"hansen");
        e_ar1_betas(r,4,1) = bic(u_hat,X,"damian");
        % Getting SE from HR and SCPC estimatots
        se_hr = HR_var(u_hat,X,X);
%         [se_scpc, cv_scpc] = scpc_var(u_hat,beta_hat,s,X,rho_bar,0.95,0);
        scpc_res = cscpc(y,X,s,0,rho_bar,0.95);
        se_scpc = scpc_res.se_beta_hat_scpc;
        cv_scpc = scpc_res.cv_scpc;
        cv_cscpc = scpc_res.cv_cscpc;
        % Testing HR and SCPC
        if m==2
            H_0 = (beta_hat - [0; beta(2:k)]);
        else
            H_0 = (beta_hat - beta);
        end
        t_stat_hr = H_0./se_hr;
        rej_hr = abs(t_stat_hr) > 1.96;
        rej_freq_m_i(r,:,1) = rej_hr;
        
        t_stat_array(r,:,1) = t_stat_hr;
        num_array(r,:,1) = H_0;
        den_array(r,:,1) = se_hr;

        t_stat_scpc = H_0./se_scpc;
        rej_scpc = abs(t_stat_scpc) > abs(cv_scpc);
        rej_freq_m_i(r,:,2) = rej_scpc;
        
        t_stat_array(r,:,2) = t_stat_scpc;
        num_array(r,:,2) = H_0;
        den_array(r,:,2) = se_scpc;

        rej_cscpc = abs(t_stat_scpc) > abs(cv_cscpc);
        rej_freq_m_i(r,:,3) = rej_cscpc;

        % Kernel for l\in L
        i=4;
        for l=L
            se_kernel = kernel_var(u_hat,X,X,D_mat,l,cholesky_flag);
            t_stat_kernel = H_0./se_kernel;
            rej_kernel = abs(t_stat_kernel) > 1.96;
            rej_freq_m_i(r,:,i) = rej_kernel;
            
            t_stat_array(r,:,i) = t_stat_kernel;
            num_array(r,:,i) = H_0;
            den_array(r,:,i) = se_kernel;
            i=i+1;
        end
        % Kernel + BSplines
        j=i; 
        n=2;
%         if r==1
%             cutoff_dist_vec = [];
%         end
        for o=1:N_order_splines
            %find optimal number of splines using Hansen's BIC
            if grid==0
                [x, fval] = fminbnd(@(x)get_opt_splines(x,s,o,y,X,cholesky_flag),4,90);
                splines_array(r,o) = round(x);
            else            
                [x, fval] = get_opt_splines_grid(s,o,y,X,cholesky_flag);        
                splines_array(r,o) = x;
            end

            for q=splines_array(r,o)%R %Number of splines

                % Matrix of BSplines on locations
                [S, delta, n_dropped]=get_bspline_mat(s,q,o,1);
                e_ar1_betas(r,5,n)=n_dropped;
                splines_array(r,o+2) = 2*delta*o;
%                 if r==1
%                     cutoff_dist_vec=[cutoff_dist_vec; NaN; L'; 2*delta*o];
%                 end
                % OLS regressors + BSplines, NO Constant
                X_bs = [X(:,2:end) S];
%                 X_bs = [X S];
                % OLS
                if m==2
                    [beta_hat_bs, u_hat_bs] = ols(y-mean(y),X_bs,X_bs, cholesky_flag); 
                else
                    [beta_hat_bs, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag);
                end
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
                e_ar1_betas(r,3,n) = fval;%bic(u_hat_bs,X_bs,"hansen");
                e_ar1_betas(r,4,n) = bic(u_hat_bs,X_bs,"damian");
                j=j+1;
                for p=[L splines_array(r,o+2)]%Cutoff lengths for Kernel
%                     disp(q);
                    % Kernel S.E. for different length cutoffs
                    se_kernel_bs = kernel_var(u_hat_bs,X_bs,X_bs,D_mat,p,cholesky_flag);
                    % t-Stat
                    % Exclude spline hat betas from estimation. 
                    % Exclude intercept from true beta
                    H_0_bs = (beta_hat_bs(1:k-1)-beta(2:k));
                    t_stat_kernel_bs = H_0_bs./se_kernel_bs(1:k-1);
                    % Rejection rule
                    rej_kernel_bs = abs(t_stat_kernel_bs) > 1.96;
                    % Collecting decisions
                    rej_freq_m_i(r,:,j) = [NaN rej_kernel_bs];
                                     
                    t_stat_array(r,:,j) = [NaN t_stat_kernel_bs];
                    num_array(r,:,j) = [NaN H_0_bs];
                    den_array(r,:,j) = [NaN se_kernel_bs(1:k-1)];
                    
                    j=j+1; %advancing rej freq index 
                end
                % advancing error AR(1) beta matrix index
                n=n+1;
            end
        end
    end
end

rej_freq(:,:) = sum(rej_freq_m_i,1)./N_reps;
avg_e_ar1(:,:) = sum(e_ar1_betas,1)./N_reps;
toc

%% Saving results %%
fprintf('Saving results\n');  

% load crank_bs_ear1_diff.mat

results_table = array2table(...
    rej_freq(:,:)',...
    'VariableNames', {'Cons.' 'Beta'});
results_table
ear1_table = array2table(...
    avg_e_ar1(:,:)',...
    'VariableNames', {'Cons.' 'Gamma' 'BIC_h' 'BIC_d' 'dropped'});
ear1_table
splines_table = array2table(splines_array, 'VariableNames', {'Qty. O1' 'Qty. O2' 'DD O1' 'DD O2'});

save(['bic_spline_' excercise '.mat'])

writetable(results_table,['../Products/bic_spline_sims_res_' excercise '.csv']);
writetable(ear1_table,['../Products/bic_spline_sims_ar1_' excercise '.csv']);





