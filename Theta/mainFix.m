if batchStartupOptionUsed
    addpath(genpath('./functions'))
    % addpath(genpath('./21200057'))
end

rng(549)
exercise = '8x8_morgans_locs_1000'; %'8x8_fix'
description = 'Quad Splines; No intercept; Morgans Locs; 1000 reps';
save_results = true;
plot_res = false;

%% Setting up parameters %%%%
n_reps = 10;
T = 500;
n_locations = 2;
theta = sqrt(2)/10;
rho = [0 0.2 0.4 0.6 1];%0.0:0.1:1.0;
l_cutoffs = 0.05:0.05:0.3;
beta = 0; %[0 0];
kernel = 'triangle';
splines_order = 3; % Order 1, step functions; 2, triangles; 3, quadratic
n_splines = 8;

%% Generate locations (fixed locations)

% s = rand(T,n_locations);
s = readmatrix("R-Morgan/coords.csv"); % Using Morgan's coordinates


%% Monte Carlo

raw_data_k_u=NaN(n_reps,length(rho),length(l_cutoffs),length(beta));
raw_data_bs_u=NaN(n_reps,length(rho),length(l_cutoffs),length(beta));
rej_freq_k_u=NaN(length(rho),length(l_cutoffs),length(beta));
rej_freq_bs_u=NaN(length(rho),length(l_cutoffs),length(beta));
n_imaginary_u=zeros(n_reps,length(rho),length(l_cutoffs),length(beta));
img_freq_bs_u=zeros(length(rho),length(l_cutoffs),length(beta));

raw_data_k=NaN(n_reps,length(rho),length(l_cutoffs),length(beta));
raw_data_bs=NaN(n_reps,length(rho),length(l_cutoffs),length(beta));
rej_freq_k=NaN(length(rho),length(l_cutoffs),length(beta));
rej_freq_bs=NaN(length(rho),length(l_cutoffs),length(beta));
n_imaginary=zeros(n_reps,length(rho),length(l_cutoffs),length(beta));
img_freq_bs=zeros(length(rho),length(l_cutoffs),length(beta));

sims_stats_mat=NaN(n_reps,length(rho),6,2);
sims_stats=NaN(length(rho),6,2);

% CI
raw_data_k_u_ci=NaN(n_reps,length(rho),length(l_cutoffs),length(beta));
raw_data_bs_u_ci=NaN(n_reps,length(rho),length(l_cutoffs),length(beta));
rej_freq_k_u_ci=NaN(length(rho),length(l_cutoffs),length(beta));
rej_freq_bs_u_ci=NaN(length(rho),length(l_cutoffs),length(beta));


raw_data_k_ci=NaN(n_reps,length(rho),length(l_cutoffs),length(beta));
raw_data_bs_ci=NaN(n_reps,length(rho),length(l_cutoffs),length(beta));
rej_freq_k_ci=NaN(length(rho),length(l_cutoffs),length(beta));
rej_freq_bs_ci=NaN(length(rho),length(l_cutoffs),length(beta));


tic

%% loop 

fprintf('Running simulations\nExcercise %s\nDescription: %s\n',exercise, description);

for r=1:n_reps

    if n_reps > 100 && mod(r,100) == 0
        fprintf('rep %d \n',r);
    elseif n_reps <= 100 && mod(r,10)==0 
        fprintf('rep %d \n',r);
    end

    % j=1;
    % for q=rho
    for j=1:length(rho)

        % Generate data for a rho
        [y, X, h] = DGP(theta,s,rho(j));

        %% OLS - 
        [beta_hat, u_hat] = ols(y,X,X,'chol');
        sims_stats_mat(r,j,1,1) = get_nn_corr(u_hat,h,1);
        sims_stats_mat(r,j,4,1) = bic(u_hat,X,'hansen');

        % Null hypothesis
        H_0 = beta_hat - beta';

        % Testing with HR
        se_hr = HR_var(u_hat,X,X);
        
        % t_stat = H_0./se_hr;
        % rej_rule_hr = abs(t_stat) > 1.96;
        % sims_stats_mat(r,j,3,1) = rej_rule_hr;

        % Ignoring imaginary se
        if (se_hr==real(se_hr))
            sims_stats_mat(r,j,6,1) = se_hr;
            t_stat_hr = H_0/se_hr;
            rej_rule_hr = abs(t_stat_hr) > 1.96;
            sims_stats_mat(r,j,3,1) = rej_rule_hr;
        end


        %% Kernel HAC with B-Splines %%
        % Find optimal number of splines, Step functions
        % [n_splines, fval] = get_opt_splines_grid(s,splines_order,y,X,'nn','chol');
        [S, ~, n_dropped] = get_bsplines(s,n_splines,splines_order);
        sims_stats_mat(r,j,2,2) = n_dropped;
        sims_stats_mat(r,j,5,2) = n_splines;
        % X_bs
        X_bs = [X S]; % Removing intercept in GDP
        % X_bs = [X(:,2:end) S]; % Removing intercept
        % X_bs = [X S(:,1:end-1)]; % Removing one column of the splines 
                                 % to avoid colinearity with intercept 
                                 % in OLS
        %% OLS with B-Splines
        [beta_hat_bs, u_hat_bs] = ols(y,X_bs,X_bs,'chol');
        sims_stats_mat(r,j,1,2) = get_nn_corr(u_hat_bs,h,1);
        sims_stats_mat(r,j,4,2) = bic(u_hat_bs,X_bs,'hansen');

        % Null Hypothesis for OLS with B-Splines
        H_0_bs = beta_hat_bs(1:length(beta)) - beta';

        % Testing with HR
        se_hr_bs = HR_var(u_hat_bs,X_bs,X_bs);
        % t_stat_hr_bs = H_0./se_hr_bs(1:length(beta));
        % rej_rule_hr_bs = abs(t_stat_hr_bs) > 1.96;
        % sims_stats_mat(r,j,3,2) = rej_rule_hr_bs;

        % Ignoring imaginary se
        if (se_hr_bs(1)==real(se_hr_bs(1))) % No Intercept, only the slope
            sims_stats_mat(r,j,6,2) = se_hr_bs(1);
            t_stat_hr_bs = H_0_bs./se_hr_bs(1);
            rej_rule_hr_bs = abs(t_stat_hr_bs) > 1.96;
            sims_stats_mat(r,j,3,2) = rej_rule_hr_bs;
        end

        % i=1;
        % for l=l_cutoffs
        for i=1:length(l_cutoffs)
            
            %% t-Tests %%
            %%%% Kernel HAC %%%%
            % SEs - varies by cutoff
            %% Uniform Kernel HAC
            se_u = kernel_var(u_hat,X,X,h,l_cutoffs(i),s,'uniform','chol');
            se_bs_u = kernel_var(u_hat_bs,X_bs,X_bs,h,l_cutoffs(i),s,'uniform','chol');
            %% Triangle Kernel HAC
            se = kernel_var(u_hat,X,X,h,l_cutoffs(i),s,kernel,'chol');
            se_bs = kernel_var(u_hat_bs,X_bs,X_bs,h,l_cutoffs(i),s,kernel,'chol');

            % Uniform Kernel
            if (se_u==real(se_u))
                t_stat_u = H_0./se_u;
                rej_rule_u = abs(t_stat_u) > 1.96;
                raw_data_k_u(r,j,i,:) = rej_rule_u;
                raw_data_k_u_ci(r,j,i,:) = 2.*abs(se_u).*1.96;
            end

            % Triangle Kernel
            if (se==real(se) )
                t_stat = H_0./se;
                rej_rule = abs(t_stat) > 1.96;
                raw_data_k(r,j,i,:) = rej_rule;
                raw_data_k_ci(r,j,i,:) =2.*abs(se).*1.96;
            end

            %%%% Kernel HAC + B-SPlines %%%%
            % t statistic and rejection for kernel+bsplines  

            % Triangle Kernels
            for k=1:length(beta)
                if (se_bs(k)==real(se_bs(k)))
                    t_stat_bs = H_0_bs(k)/se_bs(k);
                    rej_rule_bs = abs(t_stat_bs) > 1.96;
                    raw_data_bs(r,j,i,k) = rej_rule_bs;
                    raw_data_bs_ci(r,j,i,k) = 2.*abs(se_bs(k)).*1.96; 
                    % error('se are imaginary in Kernel + Bsplines');
                else
                    % raw_data_bs(r,j,i,k) = NaN;
                    n_imaginary(r,j,i,k)=1;
                end
            end

            % Uniform Kernels
            for k=1:length(beta)
                if (se_bs_u(k)==real(se_bs_u(k)))
                    t_stat_bs_u = H_0_bs(k)/se_bs_u(k);
                    rej_rule_bs_u = abs(t_stat_bs_u) > 1.96;
                    raw_data_bs_u(r,j,i,k) = rej_rule_bs_u;
                    raw_data_bs_u_ci(r,j,i,k) = 2.*abs(se_bs_u(k)).*1.96; 
                    % error('se are imaginary in Kernel + Bsplines');
                else
                    % raw_data_bs_u(r,j,i,k) = NaN;
                    n_imaginary_u(r,j,i,k)=1;
                end
            end

            % i=i+1;
        end
        % j=j+1;
    end
end

%% Summarizing simulation results

% Simulations Statistics
sims_stats(:,:,:) = mean(sims_stats_mat,1,'omitnan');

% Kernel HAC 
% Uniform Kernel
rej_freq_k_u(:,:,:) = mean(raw_data_k_u,1,'omitnan');
rej_freq_k_u_ci(:,:,:) = mean(raw_data_k_u_ci,1,'omitnan');
% Triangle Kernel
rej_freq_k(:,:,:) = mean(raw_data_k,1,'omitnan');
rej_freq_k_ci(:,:,:) = mean(raw_data_k_ci,1,'omitnan');

% Kernel HAC + BSplines
% Uniform Kernel + B-Splines
rej_freq_bs_u(:,:,:) = mean(raw_data_bs_u,1,'omitnan');
rej_freq_bs_u_ci(:,:,:) = mean(raw_data_bs_u_ci,1,'omitnan');
img_freq_bs_u(:,:,:) = mean(n_imaginary_u,1);
% Triangle Kernel + B-Splines
rej_freq_bs(:,:,:) = mean(raw_data_bs,1,'omitnan');
rej_freq_bs_ci(:,:,:) = mean(raw_data_bs_ci,1,'omitnan');
img_freq_bs(:,:,:) = mean(n_imaginary,1);

corr=rho.*exp(-(0.1/theta));
toc
%% Getting results

% Triangle Results

% rej_freq_tab_k=array2table([rho' corr' rej_freq_k(:,:,1)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs)])

% rej_freq_tab_bs=array2table([rho' corr' rej_freq_bs(:,:,1)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs)])

rej_freq_tab_k_slope=array2table([rho' corr' rej_freq_k(:,:,1) sims_stats(:,1:4,1)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])

rej_freq_tab_bs_slope=array2table([rho' corr' rej_freq_bs(:,:,1) sims_stats(:,1:4,2)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])

% img_freq_tab_bs=array2table([rho' corr' img_freq_bs(:,:,1) sims_stats(:,1:4,2)'],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])

img_freq_tab_bs_slope=array2table([rho' corr' img_freq_bs(:,:,1) sims_stats(:,1:4,2)],...
        'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])

% Uniform

% rej_freq_tab_k_u=array2table([rho' corr' rej_freq_k_u(:,:,1)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])

% rej_freq_tab_bs_u=array2table([rho' corr' rej_freq_bs_u(:,:,1)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])

rej_freq_tab_k_slope_u=array2table([rho' corr' rej_freq_k_u(:,:,1) sims_stats(:,1:4,1)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])

rej_freq_tab_bs_slope_u=array2table([rho' corr' rej_freq_bs_u(:,:,1) sims_stats(:,1:4,2)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])

% img_freq_tab_bs_u=array2table([rho' corr' img_freq_bs_u(:,:,1)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])

img_freq_tab_bs_slope_u=array2table([rho' corr' img_freq_bs_u(:,:,1) sims_stats(:,1:4,2)],...
        'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])


% CI lengths

% Triangle

rej_freq_tab_k_slope_ci=array2table([rho' corr' rej_freq_k_ci(:,:,1) sims_stats(:,1:4,1)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])

rej_freq_tab_bs_slope_ci=array2table([rho' corr' rej_freq_bs_ci(:,:,1) sims_stats(:,1:4,2)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])

% Uniform

rej_freq_tab_k_slope_u_ci=array2table([rho' corr' rej_freq_k_u_ci(:,:,1) sims_stats(:,1:4,1)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])

rej_freq_tab_bs_slope_u_ci=array2table([rho' corr' rej_freq_bs_u_ci(:,:,1) sims_stats(:,1:4,2)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC'])


%% Saving results %%

if (save_results)
    fprintf('Saving results\n');
    
    save(['outputs/theta_' exercise '.mat'])

    % Triangle
    % writetable(rej_freq_tab_k,'outputs/rej_freq_tab_kernel.csv');
    % writetable(rej_freq_tab_bs,'outputs/rej_freq_tab_kernel_bsplines.csv');
    writetable(rej_freq_tab_k_slope,['outputs/' exercise '_rej_freq_tab_kernel_slope.csv']);
    writetable(rej_freq_tab_bs_slope,['outputs/' exercise '_rej_freq_tab_kernel_bsplines_slope.csv']);
    % writetable(img_freq_tab_b[s,'output' exercise 's/img_freq_tab_kernel_bsplines.csv']);
    writetable(img_freq_tab_bs_slope,['outputs/' exercise '_img_freq_tab_kernel_bsplines_slope.csv']);

    % Uniform
    % writetable(rej_freq_tab_k_[u,'output' exercise 's/rej_freq_tab_kernel_u.csv']);
    % writetable(rej_freq_tab_bs_[u,'output' exercise 's/rej_freq_tab_kernel_bsplines_u.csv']);
    writetable(rej_freq_tab_k_slope_u,['outputs/' exercise '_rej_freq_tab_kernel_slope_u.csv']);
    writetable(rej_freq_tab_bs_slope_u,['outputs/' exercise '_rej_freq_tab_kernel_bsplines_slope_u.csv']);
    % writetable(img_freq_tab_bs_[u,'output' exercise 's/img_freq_tab_kernel_bsplines_u.csv']);
    writetable(img_freq_tab_bs_slope_u,['outputs/' exercise '_img_freq_tab_kernel_bsplines_slope_u.csv']);

    % CI 
    % Triangle
    writetable(rej_freq_tab_k_slope_ci,['outputs/' exercise '_rej_freq_tab_kernel_slope_ci.csv']);
    writetable(rej_freq_tab_bs_slope_ci,['outputs/' exercise '_rej_freq_tab_kernel_bsplines_slope_ci.csv']);

    % Uniform
    writetable(rej_freq_tab_k_slope_u_ci,['outputs/' exercise '_rej_freq_tab_kernel_slope_u_ci.csv']);
    writetable(rej_freq_tab_bs_slope_u_ci,['outputs/' exercise '_rej_freq_tab_kernel_bsplines_slope_u_ci.csv']);

end

%% Plotting %%

if plot_res

    fprintf('Plotting\n');

    load(['outputs/theta_' exercise '.mat'])

    for i=1:length(rho)
        histogram(sims_stats_mat(:,i,5,2));
        title(['Number of Splines, rho=',num2str(rho(i)) ])
        exportgraphics(gcf,strcat(['figures/hist_n_spli_rho-' num2str(rho(i)) '_' exercise '.png']))
    end

end