if batchStartupOptionUsed
    addpath(genpath('./functions'))
    addpath(genpath('./R-Morgan'))
    % addpath(genpath('./21200057'))
end

<<<<<<< HEAD
excrs = 'rej-f-rho_0.8';
=======
% excrs = 'rho_1.0_newtau';
excrs = 'rho_'+string(rho)+'_Dsqr'; % Name of the experiment
>>>>>>> a6dad9af2823db1c2441b27c6300083808af2fd5

%% Setting up parameters %%%%

rng(549) % Set random seed for reproducibility
R = 10; %MC reps
B = 1000; % Number of simulations per MC
T = 500; % Number of observations
% n_locations = 2;
b_cand = T^ ( -1/2) .* [ -3 3];%T^(-1/2).*(-10:1:10);
theta = sqrt(2)/10;
% rho = 1.0; %0.0:0.1:1; %0.0:0.1:1.0;
l_cutoffs = 0:0.025:0.15;%0.02:0.02:0.14; %;0.1; %
PC_n = 0:10:90;%20:40:100; % Number of PCs to use
splines_order = 2; % Order 1, step functions; 2, triangles; 3, quadratic
n_splines = 10;
M = 5; % Number of workers
<<<<<<< HEAD
save_results = false; % Save results
verbose = true; % Print progress
=======
save_results = true; % Save results
verbose = false; % Print progress

fprintf('Running %s with %d observations, %d MC reps, %d simulations per MC rep\n', excrs, T, R, B);
>>>>>>> a6dad9af2823db1c2441b27c6300083808af2fd5
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
% disp(size(W))

%% parallel loop 
tic;
% MC_array = zeros(length(l_cutoffs),length(PC_n)); % Array to store results
MC_array = NaN(length(l_cutoffs),length(PC_n),R);
avg_powr_array = NaN(length(l_cutoffs),length(PC_n),R); % Array to store average power
cv_array = NaN(length(l_cutoffs),length(PC_n),R); % Array
rej_array = NaN(R,1);
cil_array = NaN(R,1);

<<<<<<< HEAD
parfor (r = 1:R, M)
=======
% parfor (r = 1:R, M)
for r = 1:R
>>>>>>> a6dad9af2823db1c2441b27c6300083808af2fd5
    if mod(r, 100) == 0 || verbose
        % Print progress every 10 iterations or if verbose is true
        fprintf('MC rep %d out of %d \n', r, R);
    end

<<<<<<< HEAD
    [MC_array(:,:,r), rej_array(r), cil_array(r), avg_powr_array(:,:,r), cv_array(:,:,r)] = sim_i(theta, s, rho, h, W, l_cutoffs, PC_n, b_cand, B); % Store results of this MC rep

=======
    [MC_array(:,:,r), rej_array(r,1), cil_array(r,1), avg_powr_array(:,:,r), cv_array(:,:,r)] = sim_i(theta, s, rho, h, W, l_cutoffs, PC_n, b_cand, B); % Store results of this MC rep
    % if sum(MC_array(:,1,r), 'omitnan') > 0
    %     break;
    % end
>>>>>>> a6dad9af2823db1c2441b27c6300083808af2fd5
end
toc;

% MC_array = MC_array ./ R; % Average over MC reps
MC_array = mean(MC_array,3, "omitnan")*100; % Average over
% avg_powr_array = mean(avg_powr_array,3, "omitnan"); % Average over MC reps
% cv_array = mean(cv_array,3, "omitnan"); % Average over MC reps
MC_tbl = array2table( [l_cutoffs' MC_array],'VariableNames', ['Bwd/PC', string(PC_n) ])
cv_tbl = array2table( [l_cutoffs' mean(cv_array,3, "omitnan")], 'VariableNames', ['Bwd/PC', string(PC_n) ])
pwr_tbl = array2table( [l_cutoffs' mean(avg_powr_array,3, "omitnan")], 'VariableNames', ['Bwd/PC', string(PC_n) ])
stats_tbl = array2table( [ mean(rej_array, 'omitnan'), mean(cil_array, 'omitnan') ], 'VariableNames' ,{'Rej F', 'CI Length'}, 'RowNames', [ 'rho='+string(rho)])
<<<<<<< HEAD

=======
mean(isnan(rej_array))
>>>>>>> a6dad9af2823db1c2441b27c6300083808af2fd5

%% Save results

% if save_results
%     save(['outputs/MC_', excrs, '.mat']);
%     writetable(MC_tbl, fullfile('outputs', ['MC_', excrs, '_', num2str(T), '_', num2str(R), '_', num2str(B), '.csv']));
%     writetable(cv_tbl, fullfile('outputs', ['cv_', excrs, '_', num2str(T), '_', num2str(R), '_', num2str(B), '.csv']));
%     writetable(pwr_tbl, fullfile('outputs', ['pwr_', excrs, '_', num2str(T), '_', num2str(R), '_', num2str(B), '.csv']));
%     writetable(stats_tbl, fullfile('outputs', ['stats_', excrs, '_', num2str(T), '_', num2str(R), '_', num2str(B), '.csv']));
% end

if save_results
<<<<<<< HEAD
    save(['outputs/MC_' excrs '.mat']);
    writetable(MC_tbl, fullfile('outputs', ['MC_', excrs, '_', num2str(T), '_', num2str(R), '_', num2str(B), '.csv']));
    writetable(cv_tbl, fullfile('outputs', ['cv_', excrs, '_', num2str(T), '_', num2str(R), '_', num2str(B), '.csv']));
    writetable(pwr_tbl, fullfile('outputs', ['pwr_', excrs, '_', num2str(T), '_', num2str(R), '_', num2str(B), '.csv']));
    writetable(stats_tbl, fullfile('outputs', ['stats_', excrs, '_', num2str(T), '_', num2str(R), '_', num2str(B), '.csv']));
=======
    save(['outputs/MC_'+excrs+'.mat']);
    writetable(MC_tbl, fullfile('outputs'+['MC_'+excrs+'_'+num2str(T)+'_'+num2str(R)+'_'+num2str(B)+'.csv']));
    writetable(cv_tbl, fullfile('outputs'+['cv_'+excrs+'_'+num2str(T)+'_'+num2str(R)+'_'+num2str(B)+'.csv']));
    writetable(pwr_tbl, fullfile('outputs'+['pwr_'+excrs+'_'+num2str(T)+'_'+num2str(R)+'_'+num2str(B)+'.csv']));
    writetable(stats_tbl, fullfile('outputs'+['stats_'+excrs+'_'+num2str(T)+'_'+num2str(R)+'_'+num2str(B)+'.csv']));
>>>>>>> a6dad9af2823db1c2441b27c6300083808af2fd5
end