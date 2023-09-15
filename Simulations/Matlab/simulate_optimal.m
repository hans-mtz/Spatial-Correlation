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
excercise = 'vd_ts_gamma';%'vd_spa_gamma';
grid = 1; % 0 for fminbnd; otherwise for 5:5:90 grid
method = 'gamma'; % 'bic' for Hansen's BIC; otherwise for min residuals' AR(1) slope 
spatial = 0; % 1: S_is ~ U(0,1); 0: S ~ evenly spaced 25 by 10 in a unit dimensional square

%% Setting up parameters
global T k beta rho_bar L
N_reps = 3;
T = 250; % obs
beta = [ 1; 1.5];
k= length(beta); % number of covariates
rho_bar = 0.03; % Average pairwise correlation
L = [0.015 0.03 0.05]; % Vector of cutoff distances
M = 0.05; % 0: for Nearest Neighbour; (0,1] for distance cutoff
cholesky_flag = 'chol';
N_order_splines = 2;
N_vec_dist = 2;
N_cutoffs = length(L);
% if spatial==1
%     a=3;
% else
%     a=1;
% end
a=3;
N_estimators = (a+N_cutoffs)+(2+N_cutoffs)*N_order_splines; 
N_ols = 1+N_order_splines;

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

%% Loop
rej_freq = NaN(k,N_estimators);
avg_cv = NaN(k,2);
avg_ci = NaN(k,N_estimators);
avg_e_ar1 = NaN(5,N_ols);
t_stat_array = NaN(N_reps,k, N_estimators);
num_array = NaN(N_reps,k, N_estimators);
den_array = NaN(N_reps,k, N_estimators);
ci_array = NaN(N_reps,k, N_estimators);
cv_array = NaN(N_reps,k,2);
rej_freq_m_i = NaN(N_reps,k,N_estimators);
e_ar1_betas = NaN(N_reps,5,N_ols);
splines_array = NaN(N_reps,N_order_splines+2);

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

%         writetable(table(y, X, s), strcat('../Stata/data_splines.csv')); % Saving

    % Running OLS
    [beta_hat, u_hat] = ols(y,X,X,cholesky_flag);

    % NN corr
%     [u_t, u_t_1] = sort_u(u_hat,s);
%     e_ar1_betas(r,1:2,1) = ols(u_t, u_t_1, u_t_1, cholesky_flag);
    e_ar1_betas(r,1,1) = get_nn_corr(u_hat,D_mat,1);
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
    H_0 = (beta_hat - beta);

    t_stat_hr = H_0./se_hr;
    rej_hr = abs(t_stat_hr) > 1.96;
    rej_freq_m_i(r,:,1) = rej_hr;
    
    t_stat_array(r,:,1) = t_stat_hr;
    num_array(r,:,1) = H_0;
    den_array(r,:,1) = se_hr;
    ci_array(r,:,1) = 2.*abs(se_hr).*1.96;

%     if spatial==1
    t_stat_scpc = H_0./se_scpc;
    rej_scpc = abs(t_stat_scpc) > abs(cv_scpc);
    rej_freq_m_i(r,:,2) = rej_scpc;

    t_stat_array(r,:,2) = t_stat_scpc;
    num_array(r,:,2) = H_0;
    den_array(r,:,2) = se_scpc;
    cv_array(r,:,1) = cv_scpc;
    ci_array(r,:,2) = 2.*abs(se_scpc).*abs(cv_scpc);

    rej_cscpc = abs(t_stat_scpc) > abs(cv_cscpc);
    rej_freq_m_i(r,:,3) = rej_cscpc;

    t_stat_array(r,:,3) = t_stat_scpc;
    num_array(r,:,3) = H_0;
    den_array(r,:,3) = se_scpc;
    cv_array(r,:,2) = cv_cscpc;
    ci_array(r,:,3) = 2.*abs(se_scpc).*abs(cv_cscpc);
    i=4;
%     else
%         i=2;
%     end

    % Kernel for l\in L
    for l=L
        se_kernel = kernel_var(u_hat,X,X,D_mat,l,cholesky_flag);
        t_stat_kernel = H_0./se_kernel;
        rej_kernel = abs(t_stat_kernel) > 1.96;
        rej_freq_m_i(r,:,i) = rej_kernel;
        
        t_stat_array(r,:,i) = t_stat_kernel;
        num_array(r,:,i) = H_0;
        den_array(r,:,i) = se_kernel;
        ci_array(r,:,i) = 2.*abs(se_kernel).*1.96;
        i=i+1;
    end

    % Kernel + BSplines
    j=i; 
    n=2;
%         if r==1
%             cutoff_dist_vec = [];
%         end
    for o=1:N_order_splines
        %find optimal number of splines using Hansen's BIC or residuals' AR(1) slope      
        [x, fval] = get_opt_splines_grid(s,o,y,X,method,cholesky_flag);        
        splines_array(r,o) = x;

        for q=splines_array(r,o)%R %Number of splines

            % Matrix of BSplines on vector of locations
            [S, delta, n_dropped] = get_bsplines(s,q,o);
            e_ar1_betas(r,5,n) = n_dropped;
            splines_array(r,o+2) = 2*delta*o;
%                 if r==1
%                     cutoff_dist_vec=[cutoff_dist_vec; NaN; L'; 2*delta*o];
%                 end
            % OLS regressors + BSplines, NO Constant
            X_bs = [X(:,2:end) S];
%                 X_bs = [X S];
            % OLS
            [beta_hat_bs, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag);

            %Null hypothesis against the truth
            H_0_bs = beta_hat_bs(1:k-1)-beta(2:k);
            % Getting SE from HR 
            se_hr_bs = HR_var(u_hat_bs,X_bs,X_bs);
            t_stat_hr_bs = H_0_bs./se_hr_bs(1:k-1);
            rej_hr_bs = abs(t_stat_hr_bs) > 1.96;
            rej_freq_m_i(r,:,j) = [NaN rej_hr_bs];

            t_stat_array(r,:,j) = [NaN t_stat_hr_bs];
            num_array(r,:,j) = [NaN H_0_bs];
            den_array(r,:,j) = [NaN se_hr_bs(1:k-1)];
            ci_array(r,:,j) = [NaN 2.*abs(se_hr_bs(1:k-1)).*1.96];
            
            % NN corr
%             [u_t_bs, u_t_1_bs] = sort_u(u_hat_bs,s);
%             e_ar1_betas(r,1:2,n) = ols(u_t_bs, u_t_1_bs, u_t_1_bs, cholesky_flag);
            e_ar1_betas(r,1,n) = get_nn_corr(u_hat_bs,D_mat,1);
            e_ar1_betas(r,3,n) = bic(u_hat_bs,X_bs,"hansen");
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
avg_se(:,:) = sum(den_array,1)./N_reps;
avg_cv(:,:) = sum(cv_array,1)./N_reps;
avg_e_ar1(:,:) = sum(e_ar1_betas,1)./N_reps;
index_array = ci_array == real(ci_array);
index_table = sum(index_array,1)~=0;
ci_array(~index_array)=NaN;
avg_ci(:,:) = sum(ci_array,1,'omitnan')./sum(index_array,1);
toc

%% Preparing results %%
fprintf('Preparing results\n');  

% load crank_bs_ear1_diff.mat

results_table = array2table(...
    [rej_freq(:,:)' avg_ci(:,:)'],...
    'VariableNames', {'Cons.' 'Beta' 'CI length Cons.' 'CI length Beta'});
results_table
ear1_table = array2table(...
    avg_e_ar1(:,:)',...
    'VariableNames', {'Cons.' 'Gamma' 'BIC_h' 'BIC_d' 'dropped'});
ear1_table
splines_table = array2table(splines_array, 'VariableNames', {'Qty. O1' 'Qty. O2' 'DD O1' 'DD O2'});

%% Saving results %%
fprintf('Saving results\n');
  
save(['opt_spline_' excercise '.mat'])

writetable(results_table,['../Products/opt_spline_sims_res_' excercise '.csv']);
writetable(ear1_table,['../Products/opt_spline_sims_ar1_' excercise '.csv']);


%% Plotting %%
% clear;
fprintf('Plotting\n');
% excercise = 'grid';%'dd_2x';%'ts_gamma';%
% load(['bic_spline_' excercise '.mat'])
histogram(e_ar1_betas(:,1,1));
exportgraphics(gcf,strcat(['../Products/hist_gamma_o1_' excercise '.png']))
histogram(e_ar1_betas(:,1,2));
exportgraphics(gcf,strcat(['../Products/hist_gamma_o2_' excercise '.png']))
% Spline hist
histogram(splines_array(:,1));
exportgraphics(gcf,strcat(['../Products/hist_spline_o1_' excercise '.png']))
histogram(splines_array(:,2));
exportgraphics(gcf,strcat(['../Products/hist_spline_o2_' excercise '.png']))




