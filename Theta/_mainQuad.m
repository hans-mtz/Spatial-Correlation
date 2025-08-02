if batchStartupOptionUsed
    addpath(genpath('./functions'))
    % addpath(genpath('./21200057'))
end

rng(549)
exercise = 'finer_grid';
save_results = true;
plot_res = true;

%% Setting up parameters %%%%
n_reps = 20;
T = 50;
n_locations = 2;
theta = sqrt(2)/10;
rho = [0.0:0.1:1.0];
l_cutoffs = [0.05:0.05:0.40];
beta = [0 0];
kernel = 'triangle';
splines_order = 3;
%% Generate locations (fixed locations)

s = rand(T,n_locations);


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

sims_stats_mat=NaN(n_reps,length(rho),3,2);
sims_stats=NaN(length(rho),3,2);
tic
fprintf('Running simulations\nExcercise %s\n',exercise);


for r=1:n_reps

    if n_reps > 100 && mod(r,100) == 0
        fprintf('rep %d \n',r);
    elseif n_reps <= 100 && mod(r,10)==0 
        fprintf('rep %d \n',r);
    end

    j=1;
    for q=rho

        % Generate data for a rho
        [y, X, h] = DGP(theta,s,q);

        %% OLS - 
        [beta_hat, u_hat] = ols(y,X,X,'chol');
        sims_stats_mat(r,j,1,1) = get_nn_corr(u_hat,h,1);

        % Null hypothesis
        H_0 = beta_hat - beta';

        %% Kernel HAC with B-Splines %%
        % Find optimal number of splines, Step functions
        [n_splines, fval] = get_opt_splines_grid(s,splines_order,y,X,'nn','chol');
        [S, ~, n_dropped] = get_bsplines(s,n_splines,splines_order);
        sims_stats_mat(r,j,2,2) = n_dropped;
        sims_stats_mat(r,j,3,2) = n_splines;
        % X_bs
        X_bs = [X S(:,1:end-1)]; % Removing one column of the splines 
                                 % to avoid colinearity with intercept 
                                 % in OLS
        %% OLS with B-Splines
        [beta_hat_bs, u_hat_bs] = ols(y,X_bs,X_bs,'chol');
        sims_stats_mat(r,j,1,2) = get_nn_corr(u_hat_bs,h,1);
        % Null Hypothesis for OLS with B-Splines
        H_0_bs = beta_hat_bs(1:length(beta)) - beta';

        i=1;
        for l=l_cutoffs
            
            %% t-Tests %%
            %%%% Kernel HAC %%%%
            % SEs - varies by cutoff
            %% Uniform Kernel HAC
            se_u = kernel_var(u_hat,X,X,h,l,s,'uniform','chol');
            se_bs_u = kernel_var(u_hat_bs,X_bs,X_bs,h,l,s,'uniform','chol');
            %% Triangle Kernel HAC
            se = kernel_var(u_hat,X,X,h,l,s,kernel,'chol');
            se_bs = kernel_var(u_hat_bs,X_bs,X_bs,h,l,s,kernel,'chol');

            % Uniform Kernel
            t_stat_u = H_0./se_u;
            rej_rule_u = abs(t_stat_u) > 1.96;
            raw_data_k_u(r,j,i,:) = rej_rule_u;

            % Triangle Kernel
            t_stat = H_0./se;
            rej_rule = abs(t_stat) > 1.96;
            raw_data_k(r,j,i,:) = rej_rule;
            
            %%%% Kernel HAC + B-SPlines %%%%
            % t statistic and rejection for kernel+bsplines  

            % Triangle Kernels
            for k=1:length(beta)
                if (se_bs(k)==real(se_bs(k)))
                    t_stat_bs = H_0_bs(k)/se_bs(k);
                    rej_rule_bs = abs(t_stat_bs) > 1.96;
                    raw_data_bs(r,j,i,k) = rej_rule_bs;
                    % error('se are imaginary in Kernel + Bsplines');
                else
                    raw_data_bs(r,j,i,k) = NaN;
                    n_imaginary(r,j,i,k)=1;
                end
            end

            % Uniform Kernels
            for k=1:length(beta)
                if (se_bs_u(k)==real(se_bs_u(k)))
                    t_stat_bs_u = H_0_bs(k)/se_bs_u(k);
                    rej_rule_bs_u = abs(t_stat_bs_u) > 1.96;
                    raw_data_bs_u(r,j,i,k) = rej_rule_bs_u;
                    % error('se are imaginary in Kernel + Bsplines');
                else
                    raw_data_bs_u(r,j,i,k) = NaN;
                    n_imaginary_u(r,j,i,k)=1;
                end
            end

            i=i+1;
        end
        j=j+1;
    end
end

%% Summarizing simulation results

% Simulations Statistics
sims_stats(:,:,:) = mean(sims_stats_mat,1);
% Uniform Kernel
rej_freq_k_u(:,:,:) = mean(raw_data_k_u,1);
% Triangle Kernel
rej_freq_k(:,:,:) = mean(raw_data_k,1);

% Uniform Kernel + B-Splines
rej_freq_bs_u(:,:,:) = mean(raw_data_bs_u,1,'omitnan');
img_freq_bs_u(:,:,:) = mean(n_imaginary_u,1);
% Triangle Kernel + B-Splines
rej_freq_bs(:,:,:) = mean(raw_data_bs,1,'omitnan');
img_freq_bs(:,:,:) = mean(n_imaginary,1);

corr=rho.*exp(-(0.1/theta));
toc
%% Getting results

% Triangle Results

% rej_freq_tab_k=array2table([rho' corr' rej_freq_k(:,:,1)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs)])

% rej_freq_tab_bs=array2table([rho' corr' rej_freq_bs(:,:,1)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs)])

rej_freq_tab_k_slope=array2table([rho' corr' rej_freq_k(:,:,2) sims_stats(:,1:2,1)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop'])

rej_freq_tab_bs_slope=array2table([rho' corr' rej_freq_bs(:,:,2) sims_stats(:,1:2,2)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop'])

% img_freq_tab_bs=array2table([rho' corr' img_freq_bs(:,:,1) sims_stats(:,1:2,2)'],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop'])

img_freq_tab_bs_slope=array2table([rho' corr' img_freq_bs(:,:,2) sims_stats(:,1:2,2)],...
        'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop'])

% Uniform

% rej_freq_tab_k_u=array2table([rho' corr' rej_freq_k_u(:,:,1)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop'])

% rej_freq_tab_bs_u=array2table([rho' corr' rej_freq_bs_u(:,:,1)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop'])

rej_freq_tab_k_slope_u=array2table([rho' corr' rej_freq_k_u(:,:,2) sims_stats(:,1:2,1)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop'])

rej_freq_tab_bs_slope_u=array2table([rho' corr' rej_freq_bs_u(:,:,2) sims_stats(:,1:2,2)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop'])

% img_freq_tab_bs_u=array2table([rho' corr' img_freq_bs_u(:,:,1)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop'])

img_freq_tab_bs_slope_u=array2table([rho' corr' img_freq_bs_u(:,:,2) sims_stats(:,1:2,2)],...
        'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop'])


% Simulations stats results

% sims_stas_tab = array2table(sims_stats(:,1:2,1), 'VariableNames', {'NN' 'Dropped'})

% sims_stas_tab_bs = array2table(sims_stats(:,1:2,2), 'VariableNames', {'NN' 'Dropped'})

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

    % Stats
    % writetable(sims_stas_tab,'outputs/sims_stas_tab.csv')

end

%% Plotting %%

if plot_res

    fprintf('Plotting\n');

    load(['outputs/theta_' exercise '.mat'])

    for i=1:length(rho)
        histogram(sims_stats_mat(:,i,3,2));
        title(['Number of Splines, rho=',num2str(rho(i)) ])
        exportgraphics(gcf,strcat(['figures/hist_n_spli_rho-' num2str(rho(i)) '_' exercise '.png']))
    end

end