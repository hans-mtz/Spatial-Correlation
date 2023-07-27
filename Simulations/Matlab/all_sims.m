% Simulations - Rejection Frequencies for Kernel Var Estimator%
% Setting up parameters
global T k beta rho_bar L
T = 250;
k= 3;
beta = [ 1; 1.5; 0.5];
rho_bar = 0.03; % Average pairwise correlation
% b = 0.08;
L = [0.01 0.03 0.05 0.07];

N_models = 2; 
% Model 1 
%               e~G_exp(c_rho_bar) X~[1 N(0,1) N(0,1)]
%               Model 2 
%               e~G_exp(c_rho_bar) 
%               X~[1 G_exp(c_rho_bar) G_exp(c_rho_bar)]

N_estimators = 2+2*length(L); %HR, Kernel, Kernel + splines, SCPC
N_reps = 250;



% Fixing locations
s = rand(T,1); % vector of locations
S=get_bspline_mat(s,4,2); % Matrix of BSplines on locations

rho = [rho_bar, 0.0];
% rej_freq = zeros(N_models,N_estimators,k);
rej_freq = NaN(k,N_estimators,N_models);
% rej_freq = cell(N_estimators,N_models,k);
tic
for m=1:N_models
%     rej_freq_m = 0;
%     rej_freq_m = zeros(N_estimators,k);
    rej_freq_m = zeros(k,N_estimators);
    

    for r=1:N_reps
        % generate data      
        [y, X, D_mat] = DGP(beta,k,s,rho_bar,m);
%         writetable(table(y, X, s), strcat('../Stata/data_',num2str(m),'.csv')); % Saving
        
        % Running OLS
        [beta_hat, u_hat] = ols(y,X,X);
        [beta_hat_bs, u_hat_bs] = ols(y,[X(:,2:end) S],[X(:,2:end) S]); % OLS + BSplines, No Constant
        

        % Getting SE from HR and SCPC estimatots
        se_hr = HR_var(u_hat,X,X);
        [se_scpc, cv_scpc] = scpc_var(u_hat,beta_hat,s,X,rho_bar,0.95,0);
        
        % Testing HR and SCPC
        H_0 = (beta_hat - beta);
        t_stat_hr = H_0./se_hr;
        rej_hr = abs(t_stat_hr) > 1.96;
%         rej_freq_m(1,:) = rej_freq_m(1,:) + rej_hr'; 
        rej_freq_m(:,1) = rej_freq_m(:,1) + rej_hr;  

        t_stat_scpc = H_0./se_scpc;
        rej_scpc = abs(t_stat_scpc) > abs(cv_scpc);
%         rej_freq_m(2,:) = rej_freq_m(2,:) + rej_scpc';
        rej_freq_m(:,2) = rej_freq_m(:,2) + rej_scpc; 

        % Kernel for l\in L
        i=3;
        for l=L
            se_kernel = kernel_var(u_hat,X,X,D_mat,l);
            t_stat_kernel = H_0./se_kernel;
            rej_kernel = abs(t_stat_kernel) > 1.96;
%             rej_freq_m(i,:) = rej_freq_m(i,:) + rej_kernel';
            rej_freq_m(:,i) = rej_freq_m(:,i) + rej_kernel;
            i=i+1;
        end
        % Kernel + BSplines
        j=2+length(L)+1;
        for p=L
            se_kernel_bs = kernel_var(u_hat_bs,X,X,D_mat,p);
            t_stat_kernel_bs = (beta_hat_bs(1:k-1)-beta(2:k))./se_kernel_bs(2:k);
            rej_kernel_bs = abs(t_stat_kernel_bs) > 1.96;
            % rej_freq_m(j,:) = rej_freq_m(j,:) + rej_kernel_bs';
            rej_freq_m(2:k,j) = rej_freq_m(2:k,j) + rej_kernel_bs;
            j=j+1;
        end
        
    end
    % rej_freq(m,:,:) = rej_freq_m./N_reps;
    rej_freq(:,:,m) = rej_freq_m./N_reps;
end
toc
rej_freq
    
results_table = array2table(...
    cat(1,rej_freq(:,:,1),rej_freq(:,:,2)));

writetable(results_table,'../Products/kernel_sims.csv');
% writetable(table(y, X, s), '../Products/example_matlab.csv'); % Saving
        