if batchStartupOptionUsed
    addpath(genpath('./functions'))
    % addpath(genpath('./21200057'))
end

rng(549)
exercise = 'Opt_PC_NN_0'; %'8x8_fix'
% Preferred specification: Triangle B Splines; Gauss Kernel; Choosing # of PCs by NN: OLS with no intercept
description = '8x8 Bsplines; Obj NN set to 0; Choosing optimal PCs of Splines by NN; Gauss Kernel; Triangle Splines; Morgans Locs; OLS NO Intercept; 100 reps (par)';
save_results = true;
plot_res = true;
fix = false; %vcov fix if negative elements in the diagonal a la Cameron, Gelbech and Miller (2011)

%% Setting up parameters %%%%
n_reps = 1000;
T = 500;
n_locations = 2;
theta = sqrt(2)/10;
rho = 0.0:0.1:1; %0.0:0.1:1.0;
l_cutoffs = 0.05:0.05:0.2;
beta = 0; %[0 0]; %Intercept or no Intercept
splines_order = 2; % Order 1, step functions; 2, triangles; 3, quadratic
n_splines = 8;
M = 4; % Number of workers
%% Generate locations (fixed locations)

% If splines are fixed and locations are fixed, they need not to be estimated everytime

% Morgan's Locations
s = readmatrix("R-Morgan/coords.csv"); % Using Morgan's coordinates

% Morgan's Splines

% S = readmatrix('morgans_splines_pc'); % 
% n_dropped = 0;

% Matlab Locations
% s = rand(T,n_locations);

% Matlab Splines
[S, ~, n_dropped] = get_bsplines(s,n_splines,splines_order);

% Fixed locations means the distance matrix does not change
h=getdistmat(s,false);

% Get PCs | W is the matrix of PCs
[~,W,~] = pca(S, 'Centered', 'off');
disp(size(W))

%% Monte Carlo

raw_data_k_u=NaN(n_reps,length(rho),length(l_cutoffs),length(beta));
raw_data_bs_u=NaN(n_reps,length(rho),length(l_cutoffs),length(beta));
n_imaginary_u=zeros(n_reps,length(rho),length(l_cutoffs),length(beta));
rej_freq_k_u=NaN(length(rho),length(l_cutoffs),length(beta));
rej_freq_bs_u=NaN(length(rho),length(l_cutoffs),length(beta));
img_freq_bs_u=zeros(length(rho),length(l_cutoffs),length(beta));

raw_data_k=NaN(n_reps,length(rho),length(l_cutoffs),length(beta));
raw_data_bs=NaN(n_reps,length(rho),length(l_cutoffs),length(beta));
n_imaginary=zeros(n_reps,length(rho),length(l_cutoffs),length(beta));
rej_freq_k=NaN(length(rho),length(l_cutoffs),length(beta));
rej_freq_bs=NaN(length(rho),length(l_cutoffs),length(beta));
img_freq_bs=zeros(length(rho),length(l_cutoffs),length(beta));

% 1= NN; 2=N Drop; 3=HR rej freq; 4=BIC ; 5=N PCs; 6= HR SE
% 1=No Splines; 2=With Splines
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
parpool(M)
parfor (r=1:n_reps,M)
       
    tsk=getCurrentTask();
    
    if n_reps > 100 && mod(r,100) == 0
        fprintf('rep %d by worker %d\n',r,tsk.ID);

    elseif n_reps <= 100 && mod(r,10)==0 
        fprintf('rep %d by worker %d\n',r,tsk.ID);
 
    end    

    [ sims_stats_mat(r,:,:,:), raw_data_k_u(r,:,:,:), raw_data_k_u_ci(r,:,:,:), raw_data_k(r,:,:,:), raw_data_k_ci(r,:,:,:), raw_data_bs(r,:,:,:), raw_data_bs_ci(r,:,:,:), raw_data_bs_u(r,:,:,:), raw_data_bs_u_ci(r,:,:,:), n_imaginary(r,:,:,:), n_imaginary_u(r,:,:,:) ] = wrapper_opt_PC(rho,l_cutoffs,beta,theta,s,W,n_dropped,fix,h);
    % [sims_stats_mat(r,:,:,:),raw_data_k_u(r,:,:,:),raw_data_k_u_ci(r,:,:,:),raw_data_k(r,:,:,:),raw_data_k_ci(r,:,:,:),raw_data_bs(r,:,:,:),raw_data_bs_ci(r,:,:,:),raw_data_bs_u(r,:,:,:),raw_data_bs_u_ci(r,:,:,:),n_imaginary(r,:,:,:),n_imaginary_u(r,:,:,:)] = wrapper(rho,l_cutoffs,beta,theta,s,splines_order);

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
HR_bs_ci=mean(2.*abs(sims_stats_mat(:,:,6,1)).*1.96,1);
corr=rho.*exp(-(0.1/theta));
toc
%% Getting results

% Triangle Results

% rej_freq_tab_k=array2table([rho' corr' rej_freq_k(:,:,1)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs)])

% rej_freq_tab_bs=array2table([rho' corr' rej_freq_bs(:,:,1)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs)])

rej_freq_tab_k_slope=array2table([rho' corr' rej_freq_k(:,:,length(beta)) sims_stats(:,1:5,1)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])

rej_freq_tab_bs_slope=array2table([rho' corr' rej_freq_bs(:,:,length(beta)) sims_stats(:,1:5,2)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])

% img_freq_tab_bs=array2table([rho' corr' img_freq_bs(:,:,1) sims_stats(:,1:5,2)'],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])

rej_freq_tab_slope = array2table([ rho' corr' rej_freq_k(:,:,length(beta)) rej_freq_bs(:,:,length(beta))], ...
    'VariableNames', ['rho' 'corr' string(l_cutoffs) string(l_cutoffs)+'bs' ])

img_freq_tab_bs_slope=array2table([rho' corr' img_freq_bs(:,:,length(beta)) sims_stats(:,1:5,2)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])

rej_freq_tab_slope_tim = array2table([ rho' corr' rej_freq_k(:,1:3,length(beta)) sims_stats(:,3,1) rej_freq_bs(:,1:3,length(beta)) sims_stats(:,3,2)], ...
    'VariableNames', ['rho' 'corr' string(l_cutoffs(1:3)) 'HR' string(l_cutoffs(1:3))+'bs' 'HRbs'] )

rej_freq_tab_slope_sims_report = array2table([ rho' corr' rej_freq_k(:,1:3,length(beta)) sims_stats(:,3,1) rej_freq_bs(:,1:3,length(beta)) sims_stats(:,[3 1 5],2)], ...
        'VariableNames', ['rho' 'corr' string(l_cutoffs(1:3)) 'HR' string(l_cutoffs(1:3))+'bs' 'HRbs' 'NN' 'PCs'] )
    
% Uniform

% rej_freq_tab_k_u=array2table([rho' corr' rej_freq_k_u(:,:,1) sims_stats(:,1:5,1)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])

% rej_freq_tab_bs_u=array2table([rho' corr' rej_freq_bs_u(:,:,1) sims_stats(:,1:5,2)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])

rej_freq_tab_k_slope_u=array2table([rho' corr' rej_freq_k_u(:,:,length(beta)) sims_stats(:,1:5,1)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])

rej_freq_tab_bs_slope_u=array2table([rho' corr' rej_freq_bs_u(:,:,length(beta)) sims_stats(:,1:5,2)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])

% img_freq_tab_bs_u=array2table([rho' corr' img_freq_bs_u(:,:,1) sims_stats(:,1:5,2)],...
%     'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])

rej_freq_tab_slope_u = array2table([ rho' corr' rej_freq_k_u(:,:,length(beta)) rej_freq_bs_u(:,:,length(beta))], ...
    'VariableNames', ['rho' 'corr' string(l_cutoffs) string(l_cutoffs)+'bs' ])

img_freq_tab_bs_slope_u=array2table([rho' corr' img_freq_bs_u(:,:,length(beta)) sims_stats(:,1:5,2)],...
        'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])


% CI lengths

% Triangle

rej_freq_tab_k_slope_ci=array2table([rho' corr' rej_freq_k_ci(:,:,length(beta)) sims_stats(:,1:5,1)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])

rej_freq_tab_bs_slope_ci=array2table([rho' corr' rej_freq_bs_ci(:,:,length(beta)) sims_stats(:,1:5,2)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])

rej_freq_tab_bs_slope_ci_tim=array2table([rho' corr' rej_freq_bs_ci(:,1:3,length(beta))  HR_bs_ci'],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs(1:3)) 'HR'] )

    
% Uniform

rej_freq_tab_k_slope_u_ci=array2table([rho' corr' rej_freq_k_u_ci(:,:,length(beta)) sims_stats(:,1:5,1)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])

rej_freq_tab_bs_slope_u_ci=array2table([rho' corr' rej_freq_bs_u_ci(:,:,length(beta)) sims_stats(:,1:5,2)],...
    'VariableNames', ['$\\rho$' 'corr' string(l_cutoffs) 'NN' 'Drop' 'HR' 'BIC' 'PCs'])


%% Saving results %%

if (save_results)
    fprintf('Saving results\n');
    
    save(['outputs/theta_' exercise '.mat'])
    % load(['outputs/theta_' exercise '.mat'])

    % Triangle
    % writetable(rej_freq_tab_k,'outputs/rej_freq_tab_kernel.csv');
    % writetable(rej_freq_tab_bs,'outputs/rej_freq_tab_kernel_bsplines.csv');
    % writetable(rej_freq_tab_k_slope,['outputs/' exercise '_rej_freq_tab_kernel_slope.csv']);
    % writetable(rej_freq_tab_bs_slope,['outputs/' exercise '_rej_freq_tab_kernel_bsplines_slope.csv']);
    % writetable(img_freq_tab_bs,['outputs/' exercise '_img_freq_tab_kernel_bsplines.csv']);
    writetable(img_freq_tab_bs_slope,['outputs/' exercise '_img_freq_tab_kernel_bsplines_slope.csv']);
    writetable(rej_freq_tab_slope,['outputs/' exercise '_rej_freq_tab_slope.csv']);
    writetable(rej_freq_tab_slope_tim,['outputs/' exercise '_rej_freq_tab_slope_tim.csv']);
    writetable(rej_freq_tab_slope_sims_report,['outputs/' exercise '_rej_freq_tab_slope_sims_report.csv']);
    % Uniform
    % writetable(rej_freq_tab_k_u,['outputs/' exercise '_rej_freq_tab_kernel_u.csv']);
    % writetable(rej_freq_tab_bs_u,['outputs/' exercise '_rej_freq_tab_kernel_bsplines_u.csv']);
    % writetable(rej_freq_tab_k_slope_u,['outputs/' exercise '_rej_freq_tab_kernel_slope_u.csv']);
    % writetable(rej_freq_tab_bs_slope_u,['outputs/' exercise '_rej_freq_tab_kernel_bsplines_slope_u.csv']);
    % writetable(img_freq_tab_bs_u,['outputs/' exercise '_img_freq_tab_kernel_bsplines_u.csv']);
    writetable(img_freq_tab_bs_slope_u,['outputs/' exercise '_img_freq_tab_kernel_bsplines_slope_u.csv']);
    writetable(rej_freq_tab_slope_u,['outputs/' exercise '_rej_freq_tab_slope_u.csv']);

    % CI 
    % Triangle
    writetable(rej_freq_tab_k_slope_ci,['outputs/' exercise '_rej_freq_tab_kernel_slope_ci.csv']);
    writetable(rej_freq_tab_bs_slope_ci,['outputs/' exercise '_rej_freq_tab_kernel_bsplines_slope_ci.csv']);
    writetable(rej_freq_tab_bs_slope_ci_tim,['outputs/' exercise '_rej_freq_tab_kernel_bsplines_slope_ci_tim.csv']);


    % Uniform
    writetable(rej_freq_tab_k_slope_u_ci,['outputs/' exercise '_rej_freq_tab_kernel_slope_u_ci.csv']);
    writetable(rej_freq_tab_bs_slope_u_ci,['outputs/' exercise '_rej_freq_tab_kernel_bsplines_slope_u_ci.csv']);

end

%% Plotting %%

if plot_res

    fprintf('Plotting\n');

    % load(['outputs/theta_' exercise '.mat'])

    for i=[3,9] %i=1:length(rho)
        % histogram(sims_stats_mat(:,i,5,2));
        % title(['Number of Splines, rho=',num2str(rho(i)) ])
        % exportgraphics(gcf,strcat(['figures/hist_n_spli_rho-' num2str(rho(i)) '_' exercise '.png']))

        histogram(sims_stats_mat(:,i,1,2),20);
        title(['NN, 8x8, rho=' num2str(rho(i))])
        exportgraphics(gcf,strcat(['figures/hist_nn_rho-' num2str(rho(i)) '_' exercise '.png']))

        
    end



end