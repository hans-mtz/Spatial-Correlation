%% Simulations - Rejection Frequencies for Kernel Var Estimator%
clear;

if batchStartupOptionUsed
    addpath(genpath('./SCPC_Matlab_Implementation_20211123'))
    addpath(genpath('./21200057'))
end
%--------------GDP-------------------%
% y=beta*X+e

% Spatial Model 
% e~G_exp(c_rho_bar) 
% X~[1 G_exp(c_rho_bar)]
% 
% Time series Model
% Model 3 Spatial - Demeaned outcome
% e_s=gamma*e_s-1+u
% X_s=gamma*X_s-1+v
% X~[1 X_s],u and v are iid ~N(0,1)

%--------------Estimators---------------%
%HR, Kernel, Kernel + splines, SCPC
rng(333) % Setting seed
excercise = 'ts';%'spa'
spatial = 0; % 1: S_is ~ U(0,1); 0: S ~ evenly spaced 25 by 10 in a unit dimensional square

%% Setting up parameters
global T k beta rho_bar L
N_reps = 1000;
T = 250; % obs
beta = [ 1; 1.5];
k= length(beta); % number of covariates
rho_bar = 0.03; % Average pairwise correlation
L = [0.015 0.03 0.05]; % Vector of cutoff distances
R = 4:2:12;%[4 12 24 36 48 64 88 102 124];
M = 0; % 0: for Nearest Neighbour; (0,1] for distance cutoff

cholesky_flag = 'chol';
N_order_splines = 2;
N_vec_dist = 2;
N_cutoffs = length(L);
N_splines = length(R);
N_estimators = (3+N_cutoffs)+(2+N_cutoffs)*N_splines*N_order_splines; 
N_ols = 1+N_splines*N_order_splines;

%% Setting locations -------
if spatial==1
    % Fixing locations
    s = rand(T,N_vec_dist); % vector of locations
else % Time series
    s_temp = [repmat((1:10:T)'./T,10,1) repelem((1:25:T)'./T,25)];
    % Jittering
    v = -0.001 + (0.001+0.001).*rand(T,2);
    s = s_temp +v;
end

rej_freq = NaN(k,N_estimators);
avg_ci = NaN(k,N_estimators);
avg_e_ar1 = NaN(4,N_ols);
ci_array = NaN(N_reps,k, N_estimators);
rej_freq_m_i = NaN(N_reps,k,N_estimators);
e_ar1_betas = NaN(N_reps,4,N_ols);

tic
fprintf('Running simulations\nExercise %s\n',excercise);

for r=1:N_reps

    if N_reps<100 || mod(r,100)==0
        fprintf('rep %d \n',r);
    end

    % generate data      
    if spatial==1  
        [y, X, D_mat] = DGP(beta,s,rho_bar,1);
    else
        D_mat = getdistmat(s,false);
        [y, X, ~] = DGP_evenly(beta, 0.79,D_mat,M);
    end
%      writetable(table(y, X, s), strcat('../Products/data_ts.csv')); % Saving
    
    % Running OLS
    [beta_hat, u_hat] = ols(y,X,X,cholesky_flag);
    

    % Residuals' AR(1)
    e_ar1_betas(r,1,1) = get_nn_corr(u_hat,D_mat,1);
    e_ar1_betas(r,2,1) = bic(u_hat,X,"hansen");
    e_ar1_betas(r,3,1) = bic(u_hat,X,"damian");
    % Getting SE from HR and SCPC estimatots
    se_hr = HR_var(u_hat,X,X);
    scpc_res = cscpc(y,X,s,0,rho_bar,0.95);
    se_scpc = scpc_res.se_beta_hat_scpc;
    cv_scpc = scpc_res.cv_scpc;
    cv_cscpc = scpc_res.cv_cscpc;
    % Testing HR and SCPC
    H_0 = (beta_hat - beta);
    
    t_stat_hr = H_0./se_hr;
    rej_hr = abs(t_stat_hr) > 1.96;
    rej_freq_m_i(r,:,1) = rej_hr;
    
    ci_array(r,:,1) = 2.*abs(se_hr).*1.96;
    t_stat_scpc = H_0./se_scpc;
    rej_scpc = abs(t_stat_scpc) > abs(cv_scpc);
    rej_freq_m_i(r,:,2) = rej_scpc;
    ci_array(r,:,2) = 2.*abs(se_scpc).*abs(cv_scpc);

    rej_cscpc = abs(t_stat_scpc) > abs(cv_cscpc);
    rej_freq_m_i(r,:,3) = rej_cscpc;
    ci_array(r,:,3) = 2.*abs(se_scpc).*abs(cv_cscpc);
    i=4;

    % Kernel for l\in L
    for l=L
        se_kernel = kernel_var(u_hat,X,X,D_mat,l,cholesky_flag);
        t_stat_kernel = H_0./se_kernel;
        rej_kernel = abs(t_stat_kernel) > 1.96;
        rej_freq_m_i(r,:,i) = rej_kernel;
        ci_array(r,:,i) = 2.*abs(se_kernel).*1.96;
        i=i+1;
    end
    % Kernel + BSplines
    j=i; %2+length(L)+1;
    n=2;
    if r==1
        cutoff_dist_vec = [];
    end
    for o=1:N_order_splines
        for q=R %Number of splines grid
             % Matrix of BSplines on vector of locations
            [S, delta, n_dropped] = get_bsplines(s,q,o); 
            e_ar1_betas(r,4,n) = n_dropped;
            if r==1
                cutoff_dist_vec=[cutoff_dist_vec; NaN; L'; 2*delta*o];
            end
            % OLS regressors + BSplines, NO Constant
            X_bs = [X(:,2:end) S]; %Taking constant off, adding splines to regressors
            % OLS
            [beta_hat_bs, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag);

            %Null hypothesis against the truth
            H_0_bs = beta_hat_bs(1:k-1)-beta(2:k);
            % Getting SE from HR
            se_hr_bs = HR_var(u_hat_bs,X_bs,X_bs);
            t_stat_hr_bs = H_0_bs./se_hr_bs(1:k-1);
            rej_hr_bs = abs(t_stat_hr_bs) > 1.96;
            rej_freq_m_i(r,:,j) = [NaN rej_hr_bs];

            ci_array(r,:,j) = [NaN 2.*abs(se_hr_bs(1:k-1)).*1.96];

            % NN corr
            e_ar1_betas(r,1,n) = get_nn_corr(u_hat_bs,D_mat,1);
            e_ar1_betas(r,2,n) = bic(u_hat_bs,X_bs,"hansen");
            e_ar1_betas(r,3,n) = bic(u_hat_bs,X_bs,"damian");
            j=j+1;
            for p=[L 2*delta*o]%Cutoff lengths for Kernel
                
                se_kernel_bs = kernel_var(u_hat_bs,X_bs,X_bs,D_mat,p,cholesky_flag);
                H_0_bs = (beta_hat_bs(1:k-1)-beta(2:k));
                t_stat_kernel_bs = H_0_bs./se_kernel_bs(1:k-1);
                rej_kernel_bs = abs(t_stat_kernel_bs) > 1.96;
                rej_freq_m_i(r,:,j) = [NaN rej_kernel_bs];
                ci_array(r,:,j) = [NaN 2.*se_kernel_bs(1:k-1).*1.96];

                j=j+1; %advancing rej freq index              
            end
            % advancing error AR(1) beta matrix index
            n=n+1;
        end
    end 
end

% Taking averages across simulations
rej_freq(:,:) = sum(rej_freq_m_i,1)./N_reps;
avg_e_ar1(:,:) = sum(e_ar1_betas,1)./N_reps;
index_array = ci_array == real(ci_array);
index_table = sum(index_array,1)~=0;
ci_array(~index_array)=NaN;
avg_ci(:,:) = sum(ci_array,1,'omitnan')./sum(index_array,1);
toc

%% Preparing results %%
fprintf('Preparing results\n');  
    
results_table = array2table(...
    [rej_freq(:,:)' avg_ci(:,:)'],...
    'VariableNames', {'Cons.' 'Beta' 'CI length Cons.' 'CI length Beta'});
results_table

ear1_table = array2table(...
    avg_e_ar1(:,:)',...
    'VariableNames', {'Gamma' 'BIC_h' 'BIC_d' 'dropped'});
ear1_table

array2table(cutoff_dist_vec)

%% Saving results %%
fprintf('Saving results\n');
  
save(['grid_spline_' excercise '.mat'])

writetable(results_table,['../Products/grid_spline_sims_res_' excercise '.csv']);
writetable(ear1_table,['../Products/grid_spline_sims_ar1_' excercise '.csv']);
writetable(array2table(cutoff_dist_vec) ,['../Products/grid_spline_sims_dist_v_' excercise '.csv']);

%% Plotting %%
% clear;
fprintf('Plotting\n');
% excercise = 'grid';%'dd_2x';%'ts_gamma';%
% load(['bic_spline_' excercise '.mat'])
histogram(e_ar1_betas(:,1,1));
exportgraphics(gcf,strcat(['../Products/hist_gamma_o1_' excercise '.png']))
histogram(e_ar1_betas(:,1,2));
exportgraphics(gcf,strcat(['../Products/hist_gamma_o2_' excercise '.png']))
