if batchStartupOptionUsed
    addpath(genpath('./functions'))
    % addpath(genpath('./21200057'))
end

excrs = 'parallel';

%% Setting up parameters %%%%
fprintf("Setting up parameters for exercice " + excrs + "\n");
rng(549) % Set random seed for reproducibility
B = 200; % Number of simulations
T = 500; % Number of observations
% n_locations = 2;
b_cand = T^ ( -1/2) .* [0 -3 3];%T^(-1/2).*(-10:1:10);
theta = sqrt(2)/10;
rho = 0.5; %0.0:0.1:1; %0.0:0.1:1.0;
l_cutoffs = 0.05:0.05:0.15;%0.02:0.02:0.14; %;0.1; %
PC_n = 40:20:100;%20:40:100; % Number of PCs to use
splines_order = 2; % Order 1, step functions; 2, triangles; 3, quadratic
n_splines = 10;
M = 6; % Number of workers
save_results = true; % Save results
verbose = false; % Print progress
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

% Generate the data
fprintf('Generating data for rho = %.2f and theta = %.3f \n', rho, theta);
[y, X, ~] = DGP(theta,s,rho,false,h);

% Find tau by QMLE
if verbose 
    fprintf('Finding tau for rho = %.2f\n', rho); 
end
% Estimating tau does not change with the number of PCs
[tau, beta_hat] = get_tau(y, X, h); % epsilon_hat = y - \beta_hat* X_bs', not epsilon_hat=y-\beta_hat X-\gama_hat W; % Get tau and beta_hat using QMLE

% Simulate epsilon
[~, eps_sim] = sim_w_tau(tau, h, X, beta_hat,B);

%% Declare arrays to store results

% cv_array = NaN(B,length(l_cutoffs), length(PC_n), 1);
tstat_array = NaN(B,length(l_cutoffs), length(PC_n), length(b_cand));
se_array = NaN(B, length(l_cutoffs), length(PC_n), length(b_cand)); % Standard errors for t-statistics

%% Run the simulations in parallel
loopStart = tic; % Start timer for the entire simulation

parpool(M); % Start parallel pool with M workers
parfor (i = 1:B, M)
    fprintf('Starting simulation %d of %d \n', i, B);

    % Call the wrapper function to perform cross-validation
    [tstat_array(i,:,:,:), se_array(i,:,:,:)] = wrapperCV(PC_n, W, X, eps_sim(:,i), beta_hat, h, s, b_cand, l_cutoffs, false);
end
loopEnd = toc(loopStart);
fprintf('Total time for simulations: %.2f seconds\n', loopEnd);

%% Reshape the results to match the expected output format
cv_array = NaN(length(l_cutoffs), length(PC_n)); % Lower and upper bounds for CV
power_array = NaN(length(l_cutoffs), length(PC_n), length(b_cand)-1);

loopCV = tic; % Start timer for critical values
% Get critical values for each pair of l_cutoff and number of PCs
for k = 1:length(b_cand)
    fprintf('Calculating critical values for beta candidate %d of %d \n', k, length(b_cand));
    % Reshape tstat_array to match the expected dimensions
    tstat_i = squeeze(tstat_array(:,:,:,k));
    % Calculate the critical values for each l_cutoff and number of PCs
    if k > 1, kk= k-1; end
    for i = 1:length(PC_n)
        for j = 1:length(l_cutoffs)
            dist_i = fitdist(abs(tstat_i(:,j,i)),"Kernel","Kernel","normal");
            % Calculate the critical values for each beta candidate
            if k==1 
                cv_array(j,i) = icdf(dist_i, 0.95);
            else
                power_array(j,i,kk) = 1 - cdf(dist_i, cv_array(j,i));
            end
        end
    end
end

loopCVEnd = toc(loopCV);
fprintf('Total time for critical values: %.2f seconds\n', loopCVEnd);
bothEnd = toc(loopStart);
fprintf('Total time : %.2f seconds\n', bothEnd);

%% Collect results in tables

power_mean = mean(power_array,3); % Mean power for each l_cutoff and number of PCs

if size(power_mean,2) <= 1

    bw_tbl = array2table( [l_cutoffs' cv_array(:,:,1) cv_array(:,:,2) power_mean],'VariableNames', {'Cutoff','CV_Lower', 'CV_Upper', 'Mean Power'})
else
    bw_tbl = array2table( [l_cutoffs' power_mean],'VariableNames', ['Bwd/PC', string(PC_n) ])
    cv_tbl = array2table( [l_cutoffs' cv_array(:,:)], 'VariableNames', ['Bwd/PC', string(PC_n) ])
    % cv_u_tbl = array2table( [l_cutoffs' cv_array(:,:,2)], 'VariableNames', ['Bwd/PC', string(PC_n) ])
    % cv_l_tbl = array2table( [l_cutoffs' cv_array(:,:,1)], 'VariableNames', ['Bwd/PC', string(PC_n) ])
end


%% Save the results

if save_results
    fprintf('Saving results\n');
    
    save(['outputs/' excrs 'bw_cv_sel.mat'])

    writetable(bw_tbl,['outputs/' excrs 'bw_cv_pwr.csv']);
    writetable(cv_tbl,['outputs/' excrs 'bw_cv.csv']);
    % writetable(cv_u_tbl,['outputs/' excrs 'bw_cv_u.csv']);
    % writetable(cv_l_tbl,['outputs/' excrs 'bw_cv_l.csv']);
end

%% Close parallel pool

if ~batchStartupOptionUsed
    fprintf('Closing parallel pool\n');
    delete(gcp('nocreate')) % Close the parallel pool
end