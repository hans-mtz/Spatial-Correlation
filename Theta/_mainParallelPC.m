if batchStartupOptionUsed
    addpath(genpath('./functions'))

end

rng(549)
exercise = '8x8_by_PC'; %'8x8_fix'
description = 'By PC; Gauss Kernel; Triangle Splines; Morgans Locs; OLS NO Intercept; 1000 reps (par)';
save_results = true;
plot_res = true;
fix = false; %vcov fix if negative elements in the diagonal a la Cameron, Gelbech and Miller (2011)

%% Setting up parameters %%%%
n_reps = 1000;
T = 500;
n_locations = 2;
theta = sqrt(2)/10;
rho = 0.8; %0.0:0.1:1.0;
l_cutoffs = 0.1;
beta = 0; %[0 0]; %Intercept or no Intercept
splines_order = 2; % Order 1, step functions; 2, triangles; 3, quadratic
n_splines = 8;
M = 6;
n_pc = n_splines^2;

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

% Get PCs
[~,W,~] = pca(S, 'Centered', 'off');
% [W,~,n_dropped] = get_PC(n_pc,s,n_pc,splines_order);
disp(size(W))
%% Monte Carlo

raw_data_bs=NaN(n_reps,size(W,2),length(beta));
n_imaginary=zeros(n_reps,size(W,2),length(beta));
rej_freq_bs=NaN(size(W,2),length(beta));
img_freq_bs=zeros(size(W,2),length(beta));

% 1= NN; 2=N Drop; 3=HR rej freq; 4=BIC ; 5=N Splines; 6= HR SE
% 1=No Splines; 2=With Splines
sims_stats_mat=NaN(n_reps,size(W,2),4);
sims_stats=NaN(size(W,2),4);

% CI

raw_data_bs_ci=NaN(n_reps,size(W,2),length(beta));
rej_freq_bs_ci=NaN(size(W,2),length(beta));


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

    [ sims_stats_mat(r,:,:),  raw_data_bs(r,:,:), raw_data_bs_ci(r,:,:),  n_imaginary(r,:,:)] = wrapper_PC(rho,l_cutoffs,beta,theta,s,W,fix,h);


end

%% Summarizing simulation results

% Simulations Statistics
sims_stats(:,:) = mean(sims_stats_mat,1,'omitnan');

% Kernel HAC + BSplines

% Triangle Kernel + B-Splines
rej_freq_bs(:,:) = mean(raw_data_bs,1,'omitnan');
rej_freq_bs_ci(:,:) = mean(raw_data_bs_ci,1,'omitnan');
img_freq_bs(:,:) = mean(n_imaginary,1);
HR_bs_ci=mean(2.*abs(sims_stats_mat(:,:,4)).*1.96,1);
corr=rho.*exp(-(0.1/theta));
toc

%% Getting results


%HAC 0.1 Rej Freq
tbl_names = {'HAC0.1', 'CI', 'HR', 'BIC', 'NN'}
pc_tbl = array2table( [rej_freq_bs rej_freq_bs_ci sims_stats(:,1:3)], 'VariableNames',tbl_names)

fig_names = {'HAC 2\sigma = .1 Reject', 'CI Length', 'HR Reject', 'BIC', 'NN Corr'};

%% Saving results %%

if (save_results)
    fprintf('Saving results\n');
    
    save(['outputs/theta_' exercise '.mat'])

    writetable(pc_tbl,['outputs/' exercise '_pc_plot_data.csv']);
    
end

%% Plotting %%

if plot_res

    fprintf('Plotting\n');

    for i = 1:size(pc_tbl, 2)
        subplot(ceil(sqrt(size(pc_tbl, 2))), ceil(sqrt(size(pc_tbl, 2))), i);
        plot(pc_tbl, tbl_names{i},'Color','black');
        ylim padded
        ylabel('')
        title(fig_names{i});
        if i==1
            ylim([0.05 0.275]);
            yticks([0.05:0.05:0.25]);
        end

        if i==size(pc_tbl, 2)
            % ylim([0.05 0.275]);
            yticks([0.1:0.2:0.5]);
        end
    end

    sgtitle('Properties of Alternate G Specifications Using Principal Components');
    xlabel('Number of Principal Components');

    exportgraphics(gcf,strcat(['figures/' exercise '_plot.pdf']))

end