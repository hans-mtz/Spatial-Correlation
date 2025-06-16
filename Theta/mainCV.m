if batchStartupOptionUsed
    addpath(genpath('./functions'))
    % addpath(genpath('./21200057'))
end

rng(549)
exercise = 'Opt_PC_NN_0'; %'8x8_fix'
% Preferred specification: Triangle B Splines; Gauss Kernel; Choosing # of PCs by NN: OLS with no intercept
description = '';
save_results = true;
plot_res = true;

%% Setting up parameters %%%%

B = 100; % Number of simulations
% T = 500; % Number of observations
% n_locations = 2;
theta = sqrt(2)/10;
rho = 0.5; %0.0:0.1:1; %0.0:0.1:1.0;
l_cutoffs = 0.1; %0.05:0.05:0.2;
splines_order = 2; % Order 1, step functions; 2, triangles; 3, quadratic
n_splines = 10;
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

% Generate the data

[y, X, h] = DGP(theta,s,rho,false,h);



%% For a given rho, l_cutoff, and number of PCs from the Tensor product of the Bsplines

PC_n = 20;
X_bs = [X W(:,1:PC_n)]; % X with B-splines, no intercept

% Find tau by QMLE
[tau, beta_hat] = get_tau(y, X_bs, h);


%% Simulate data given tau, for beta_cand

% For C.V. beta candidate is always zero

beta_cand = [ 0 beta_hat(2:end)']; % Candidate beta and estimates of the Bsplines coefficients

%% Simulate data with the candidate beta and the estimated tau

[y_sim, eps_sim] = sim_w_tau(tau, h, X_bs, beta_cand,B);

% Get the statistic for the simulated data
stat_vec = NaN(B,1); % Store the distribution of the statistic
se_vec = NaN(B,1); % Store the standard errors of the statistic
for i = 1:B
    [beta_hat_sim, eps_hat_sim] = ols(y_sim(:,i), X_bs, X_bs, 'chol');
    % For CV Test H_0: beta_0 = 0 | beta =0
    [t_val, se_kvar] = get_statistic(beta_hat_sim, zeros(length(beta_hat_sim),1), eps_hat_sim, X_bs, h, l_cutoffs, s);
    stat_vec(i) = t_val(1);
    se_vec(i) = se_kvar(1);
end

%%  Get the empirical distribution of the statistic
x = linspace(min(stat_vec), max(stat_vec), 100);
% stat_density = paretotails(stat_vec, 0.01, 0.99);
stat_density = fitdist(stat_vec, "Kernel", "Kernel", "normal");

% [f, x1] = ecdf(stat_vec);
% type1_density = makedist('PiecewiseLinear', 'x', x1, 'Fx', f); 
ecdf(stat_vec);
hold on
plot(x,cdf(stat_density,x),'r--');
% hold on
% plot(x, cdf(type1_density,x),'b-');

%% Get critival values for the statistic
cv_h = icdf(stat_density,[0.05 0.95]); % 95% quantile of the statistic distribution

%%  Get the power of the test

%for a candidate of beta

beta_cand_pwr = [ 0.5 beta_hat(2:end)']; % Candidate beta and estimates of the Bsplines coefficients

%% Simulate data with the candidate beta and the estimated tau

% B = 100; % Number of simulations
[y_sim_pwr, eps_sim_pwr] = sim_w_tau(tau, h, X_bs, beta_cand_pwr,B);

% Get the statistic for the simulated data
pwr_stat_vec = NaN(B,1); % Store the distribution of the statistic
pwr_se_vec = NaN(B,1); % Store the standard errors of the statistic
for i = 1:B
    [beta_hat_sim_cand, eps_hat_sim_cand] = ols(y_sim_pwr(:,i), X_bs, X_bs, 'chol');
    % For CV Test H_0: beta_0 = 0 | beta = beta_cand
    [pwr_t_val, pwr_se_kvar] = get_statistic(beta_hat_sim_cand, zeros(length(beta_hat_sim_cand),1), eps_hat_sim_cand, X_bs, h, l_cutoffs, s);
    pwr_stat_vec(i) = pwr_t_val(1);
    pwr_se_vec(i) = pwr_se_kvar(1);
end

%%  Get the empirical distribution of the statistic

x = linspace(min(pwr_stat_vec), max(pwr_stat_vec), 100);
% pwr_stat_density = paretotails(pwr_stat_vec, 0.01, 0.99);
pwr_stat_density = fitdist(pwr_stat_vec, "Kernel", "Kernel", "normal");

% [f, x1] = ecdf(pwr_stat_vec);
% type2_density = makedist('PiecewiseLinear', 'x', x1, 'Fx', f); 
ecdf(pwr_stat_vec);
hold on
plot(x,cdf(pwr_stat_density,x),'r--');
% hold on
% plot(x, cdf(type2_density,x),'b-');

%% Get the probability of rejection for the candidate beta
% using CV from testing the null given beta_cand =0

probs = cdf(pwr_stat_density, cv_h); % Probability of rejection for the candidate beta
prob_rejection = probs(1) + (1 - probs(2)); % Probability of rejection for the candidate beta

