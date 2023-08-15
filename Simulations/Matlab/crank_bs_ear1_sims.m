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
T = 2500; % obs
beta = [ 1; 1.5];
k= length(beta); % number of covariates
rho_bar = 0.03; % Average pairwise correlation
L = [0.015 0.03 0.05]; % Vector of cutoff distances
R = [4 12 24 36];% 90 120 180]; % Number of triangles
% R = 40:48; % 43 for order 1 splines
% R = 64:84; % 82 for order 2 splines
cholesky_flag = 'chol';
N_order_splines = 2;
N_models = 2;
N_reps = 10;

N_cutoffs = length(L);
N_triangles = length(R);
N_estimators = 2+N_cutoffs+N_cutoffs*N_triangles*N_order_splines; 
N_ols = 1+N_triangles*N_order_splines;

% Fixing locations
s = rand(T,1); % vector of locations
rej_freq = NaN(k,N_estimators,N_models);
avg_e_ar1 = NaN(2,N_ols, N_models);
t_stat_array = NaN(N_reps,k, N_estimators, N_models);
num_array = NaN(N_reps,k, N_estimators, N_models);
den_array = NaN(N_reps,k, N_estimators, N_models);

tic
% for m=2:3
    
rej_freq_m_i = NaN(N_reps,k,N_estimators);
e_ar1_betas = NaN(N_reps,2,N_ols);

    for r=1:N_reps
        
        % generate data      
        [y, X, D_mat] = DGP(beta,s,rho_bar,2);

%         writetable(table(y, X, s), strcat('../Stata/data_splines.csv')); % Saving
        for m=2:3
            fprintf('model %d, rep %d \n',m,r);
            % Running OLS
            if m==3
                [beta_hat, u_hat] = ols(y-mean(y),X,X,cholesky_flag);
            else
                [beta_hat, u_hat] = ols(y,X,X,cholesky_flag);
            end
    
            % Residuals' AR(1)
            [u_t, u_t_1] = sort_u(u_hat,s);
            e_ar1_betas(r,:,1) = ols(u_t, u_t_1, u_t_1, cholesky_flag);
            
            % Getting SE from HR and SCPC estimatots
            se_hr = HR_var(u_hat,X,X);
            [se_scpc, cv_scpc] = scpc_var(u_hat,beta_hat,s,X,rho_bar,0.95,0);
            
            % Testing HR and SCPC
            if m==3
            H_0 = (beta_hat - [0; beta(2:k)]);
        else
            H_0 = (beta_hat - beta);
            end
        t_stat_hr = H_0./se_hr;
            rej_hr = abs(t_stat_hr) > 1.96;
            rej_freq_m_i(r,:,1) = rej_hr;
            
            t_stat_array(r,:,1,m-1) = t_stat_hr;
            num_array(r,:,1,m-1) = H_0;
            den_array(r,:,1,m-1) = se_hr;
    
            t_stat_scpc = H_0./se_scpc;
            rej_scpc = abs(t_stat_scpc) > abs(cv_scpc);
            rej_freq_m_i(r,:,2) = rej_scpc;
            
            t_stat_array(r,:,2,m-1) = t_stat_scpc;
            num_array(r,:,2,m-1) = H_0;
            den_array(r,:,2,m-1) = se_scpc;
    
            % Kernel for l\in L
            i=3;
            for l=L
                se_kernel = kernel_var(u_hat,X,X,D_mat,l, cholesky_flag);
                t_stat_kernel = H_0./se_kernel;
                rej_kernel = abs(t_stat_kernel) > 1.96;
                rej_freq_m_i(r,:,i) = rej_kernel;
                
                t_stat_array(r,:,i,m-1) = t_stat_kernel;
                num_array(r,:,i,m-1) = H_0;
                den_array(r,:,i,m-1) = se_kernel;
                i=i+1;
            end
            % Kernel + BSplines
            j=i; %2+length(L)+1;
            n=2;
            for o=1:N_order_splines
                for q=R %Number of splines
    
                    % Matrix of BSplines on locations
                    S=get_bspline_mat(s,q,o);
                    % OLS regressors + BSplines, NO Constant
                    X_bs = [X(:,2:end) S];
    %                 X_bs = [X S];
                    % OLS
                    if m==3
                        [beta_hat_bs, u_hat_bs] = ols(y-mean(y),X_bs,X_bs, cholesky_flag); 
                    else
                        [beta_hat_bs, u_hat_bs] = ols(y,X_bs,X_bs, cholesky_flag);
                    end
                    % Residuals' AR(1)
                    [u_t_bs, u_t_1_bs] = sort_u(u_hat_bs,s);
                    e_ar1_betas(r,:,n) = ols(u_t_bs, u_t_1_bs, u_t_1_bs, cholesky_flag);
    
                    for p=L %Cutoff lengths for Kernel
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
                        
                        
                        t_stat_array(r,:,j,m-1) = [NaN t_stat_kernel_bs];
                        num_array(r,:,j,m-1) = [NaN H_0_bs];
                        den_array(r,:,j,m-1) = [NaN se_kernel_bs(1:k-1)];
                        
                        j=j+1;
                        
                    end
                    % advancing error AR(1) beta matrix index
                    n=n+1;
                end
            end
            rej_freq(:,:,m-1) = sum(rej_freq_m_i,1)./N_reps;
            avg_e_ar1(:,:,m-1) = sum(e_ar1_betas,1)./N_reps;
        end
    end

% end
toc

% %% Saving results %%
%     
results_table = array2table(...
    cat(2,rej_freq(:,:,1)',rej_freq(:,:,2)'));
results_table
ear1_table = array2table(...
    cat(2,avg_e_ar1(:,:,1)',avg_e_ar1(:,:,2)'));
ear1_table

save crank_bs_ear1_2500.mat

writetable(results_table,'../Products/crank_bs_ear1_sims_res_2500.csv');
writetable(ear1_table,'../Products/crank_bs_ear1_sims_ar1_2500.csv');

%% Plotting %%
% 
% % load('crank_bs_ear1.mat')
% 
% for j=1:length(t_stat_array(1,1,:,1))
% 
%     select = t_stat_array(:,2,j,1) == real(t_stat_array(:,2,j,1));
%     
%     histogram(t_stat_array(select,2,j,1));
%     exportgraphics(gcf,strcat('../Products/hist_t_',num2str(j),'.png'))
%     
%     select = num_array(:,2,j,1) == real(num_array(:,2,j,1));
%     
%     h=histogram(num_array(select,2,j,1));
%     h.FaceColor = [240/255 128/255 128/255];
%     exportgraphics(gcf,strcat('../Products/hist_num_',num2str(j),'.png'))
%     
%     select = den_array(:,2,j,1) == real(den_array(:,2,j,1));
%     
%     h=histogram(den_array(select,2,j,1));
%     h.FaceColor = [102/255 205/255 170/255];
%     exportgraphics(gcf,strcat('../Products/hist_den_',num2str(j),'.png'))
%     
% end
% 
% 
% %% comples SE
% 
% cmlpx_how=zeros(29,1);
% 
% for j=1:29
%     cmlpx_how(j) = sum(t_stat_array(:,2,j,1) ~= real(t_stat_array(:,2,j,1)));
% end
% 
% 
% select = t_stat_array(:,2,j,1) == real(t_stat_array(:,2,j,1));

