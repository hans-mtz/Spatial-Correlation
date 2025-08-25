function [id_mat, rej_rule, ci_length , power_mean, cv_array] = sim_i(theta, s, rho, h, W, l_cutoffs, PC_n, b_cand, B)
    
T = size(s,1); % Number of observations
% if batchStartupOptionUsed
%     addpath(genpath('./functions'))
%     addpath(genpath('./R-Morgan'))
%     % addpath(genpath('./21200057'))
% end

% excrs = 'Fast-with-df-1MC';

% %% Setting up parameters %%%%

% rng(549) % Set random seed for reproducibility
% B = 1000; % Number of simulations
% T = 500; % Number of observations
% % n_locations = 2;
% b_cand = T^ ( -1/2) .* [ -3 3];%T^(-1/2).*(-10:1:10);
% theta = sqrt(2)/10;
% rho = 0.5; %0.0:0.1:1; %0.0:0.1:1.0;
% l_cutoffs = 0:0.025:0.15;%0.02:0.02:0.14; %;0.1; %
% PC_n = 0:10:90;%20:40:100; % Number of PCs to use
% splines_order = 2; % Order 1, step functions; 2, triangles; 3, quadratic
% n_splines = 10;
% M = 5; % Number of workers
save_results = false; % Save results
verbose = false; % Print progress
% %% Generate locations (fixed locations)

% % If splines are fixed and locations are fixed, they need not to be estimated everytime

% % Morgan's Locations
% s = readmatrix("R-Morgan/coords.csv"); % Using Morgan's coordinates

% % Morgan's Splines

% % S = readmatrix('morgans_splines_pc'); % 
% % n_dropped = 0;

% % Matlab Locations
% % s = rand(T,n_locations);

% % Matlab Splines
% [S, ~, n_dropped] = get_bsplines(s,n_splines,splines_order);

% % Fixed locations means the distance matrix does not change
% h=getdistmat(s,false);

% % Get PCs | W is the matrix of PCs
% [~,W,~] = pca(S, 'Centered', 'off');
% % disp(size(W))

% Generate the data
if verbose 
    fprintf('Generating data for rho = %.2f and theta = %.3f \n', rho, theta); 
end
[y, X, ~] = DGP(theta,s,rho,false,h);

% Find tau by QMLE
if verbose 
    fprintf('Finding tau for rho = %.2f, and %d PCs\n', rho, ii); 
end
% Estimating tau does not change with the number of PCs
[tau, beta_hat, exitflag] = get_tau(y, X, h, false); % epsilon_hat = y - \beta_hat* X_bs', not epsilon_hat=y-\beta_hat X-\gama_hat W; % Get tau and beta_hat using QMLE

if exitflag <= 0
    id_mat = NaN(length(l_cutoffs), length(PC_n));
    power_mean = NaN(length(l_cutoffs), length(PC_n));
    cv_array = NaN(length(l_cutoffs), length(PC_n), 1);
    rej_rule = NaN;
    ci_length = NaN;
<<<<<<< HEAD
    warning('Optimization did not converge, exiting simulation');
=======
    warning('Optimization did not converge, exiting simulation. Tau values: %f, %f', tau(1), tau(2));
>>>>>>> a6dad9af2823db1c2441b27c6300083808af2fd5
    return
end

%% Declare arrays to store results

cv_array = NaN(length(l_cutoffs), length(PC_n), 1);
power_array = NaN(length(l_cutoffs), length(PC_n), length(b_cand));

%% For a given rho, l_cutoff, and number of PCs from the Tensor product of the Bsplines
j = 0;
jj = 0;
k = 0;
kk = 0;

loopStart = tic;

% % % parpool(M)

for sp = PC_n

    if verbose 
        fprintf('Starting Simulatios for %d Number of PCs \n', sp); 
    end
    k = k + 1;
    ii = min(sp, size(W,2)); % Ensure sp does not exceed the number of PCs available
    X_bs = [X W(:,1:ii)]; % X with B-splines, no intercept
    
    %% Simulate data given tau, for beta_cand
    
    % For C.V. beta, candidate is always zero
    if length(beta_hat) > 1
        beta_cand = [ 0 beta_hat(2:end)']; % Candidate beta and estimates of the Bsplines coefficients
    else 
        beta_cand = 0; % If there is no B-splines, the candidate beta is zero
    end
    %% Simulate data with the candidate beta and the estimated tau
    
    % fprintf('Simulating data with tau = %.4f, %.4f and beta_cand = %.4f with %.0f repetitions \n', tau(1), tau(2), beta_cand(1),B);
    [y_sim, eps_sim] = sim_w_tau(tau, h, X, beta_cand,B);
    XX_inv = inv(X_bs'*X_bs);
    M_x = eye(T) - X_bs*XX_inv*X_bs';
    M_w = eye(T) - X_bs(:,2:end)*inv(X_bs(:,2:end)'*X_bs(:,2:end))*X_bs(:,2:end)'; %Form Projection matrix away from W
    Xtilde = M_w*X; %Partial out W from X
    
    cvTime = tic; % Start timer for critical values
    % Get the statistic for the simulated data
    stat_vec = NaN(B,1); % Store the distribution of the statistic
    se_vec = NaN(B,1); % Store the standard errors of the statistic
    j= 0;
    for l = l_cutoffs
        j = j + 1;
        % fprintf('Finding critical values for %d PCs and %0.4f cutoff', sp, l);
        % % % parfor (i = 1:B, M)
        % % % for i = 1:B
            % [beta_hat_sim, eps_hat_sim] = ols(y_sim(:,i), X_bs, X_bs, 'chol');
            eps_hat_sim = M_x*y_sim;
            beta_hat_sim = XX_inv'*(X_bs'*y_sim);
            b_k = size(beta_hat_sim,1);
            % For CV Test H_0: beta_0 = 0 | beta =0
            S = Xtilde.*eps_hat_sim; %Form Scores
            % K = triangle_weights(s,l);  %Form weight matrixDKD
            if l <= 0.0 
                K = eye(T);
            else 
                K = normpdf(h,0,l/2)./normpdf(0,0,l/2);
            end
            se_num  = diag(S'*K*S).^.5;   %Numerator
            se_denom = sum(Xtilde.^2);  %Denominator
            se_vec = (se_num./se_denom)*sqrt(T/(T-b_k));% .* sqrt(T/(T-size(beta_hat_sim,2)));  %Standard errors
            stat_vec= (beta_hat_sim(1,:)./se_vec')'; %tstats 
            %[t_val, se_kvar] = get_statistic(beta_hat_sim, eps_hat_sim, X_bs, h, l, s);
            %stat_vec(i) = t_val(1);
            %se_vec(i) = se_kvar(1);
     

%  Get the empirical distribution of the statistic

        % stat_density = paretotails(stat_vec, 0.01, 0.99);
        % Fit a kernel density to the statistic
        % Use 95% of the absolute values of the statistic
        stat_density = fitdist(abs(stat_vec), "Kernel", "Kernel", "normal");

        % Plotting
        % x = linspace(min(stat_vec), max(stat_vec), 100);
        % % [f, x1] = ecdf(stat_vec);
        % % type1_density = makedist('PiecewiseLinear', 'x', x1, 'Fx', f); 
        % ecdf(stat_vec);
        % hold on
        % plot(x,cdf(stat_density,x),'r--');
        % hold on
        % plot(x, cdf(type1_density,x),'b-');

%% Get critival values for the statistic
        % cv_h = icdf(stat_density,[0.025 0.975]); % 95% quantile of the statistic distribution
        % fprintf('Critical values for %d PCs and %0.4f cutoff: %.4f, %.4f \n', sp, l, cv_h(1), cv_h(2));
        % cv_array(j,k,:) = cv_h; % Store the critical values for the statistic
        cv_h = max(1.96, icdf(stat_density, 0.95)); % 95% quantile of the abs statistic distribution
        if verbose 
            fprintf('Critical values for %d PCs and %0.4f cutoff: %.4f \n', sp, l, cv_h);
        end
        cv_array(j,k,1) = cv_h; % Store the critical values for the statistic
%%  Get the power of the test
    end
    cvTimeEnd = toc(cvTime);
    if verbose
        fprintf('Finished critical values for %d PCs in %.2f seconds \n', sp, cvTimeEnd);
    end
    kk = 0; % Reset kk for the next candidate beta

    pwrTime = tic; % Start timer for power simulations
    for b = b_cand
        if verbose 
            fprintf('Power: Starting Simulations for %d Number of PCs and beta_cand = %.4f \n', sp, b);
        end
        %for a candidate of beta

        if (length(beta_hat) > 1)  
            beta_cand_pwr = [ b beta_hat(2:end)']; % Candidate beta and estimates of the Bsplines coefficients
        else
            beta_cand_pwr = b;
        end
        
        % Simulate data with the candidate beta and the estimated tau

        % B = 100; % Number of simulations
        % fprintf('Simulating data with tau = %.4f, %.4f and beta_cand = %.4f for %.0f repetitions \n', tau(1), tau(2), beta_cand_pwr(1), B);
        
        % [y_sim_pwr, eps_sim_pwr] = sim_w_tau(tau, h, X, beta_cand_pwr,B); % Drawing new eps_sim
        y_sim_pwr = X*beta_cand_pwr' + eps_sim; % Fixing eps_sim

        % Get the statistic for the simulated data
        pwr_stat_vec = NaN(B,1); % Store the distribution of the statistic
        pwr_se_vec = NaN(B,1); % Store the standard errors of the statistic
        
        kk = kk + 1;
        jj = 0;
        for l = l_cutoffs
            jj = jj + 1;
            % % % parfor (i = 1:B,M)
            %for i = 1:B
                % [beta_hat_sim_cand, eps_hat_sim_cand] = ols(y_sim_pwr(:,i), X_bs, X_bs, 'chol');
                % For CV Test H_0: beta_0 = 0 | beta = beta_cand
            %    eps_hat_sim_cand = M_x*y_sim_pwr(:,i);
            %    beta_hat_sim_cand = XX_inv'*(X_bs'*y_sim_pwr(:,i));
            %    [pwr_t_val, pwr_se_kvar] = get_statistic(beta_hat_sim_cand, eps_hat_sim_cand, X_bs, h, l, s);
            %    pwr_stat_vec(i) = pwr_t_val(1);
            %    pwr_se_vec(i) = pwr_se_kvar(1);
            %end

            % [beta_hat_sim, eps_hat_sim] = ols(y_sim(:,i), X_bs, X_bs, 'chol');
            eps_hat_sim_cand = M_x*y_sim_pwr;
            beta_hat_sim_cand = XX_inv'*(X_bs'*y_sim_pwr);
            b_k = size(beta_hat_sim_cand,1);
            % For CV Test H_0: beta_0 = 0 | beta =0
            S = Xtilde.*eps_hat_sim_cand; %Form Scores
            % K = triangle_weights(s,l);  %Form weight matrix
            if l <= 0.0 
                K = eye(T);
            else 
                K = normpdf(h,0,l/2)./normpdf(0,0,l/2);
            end
            se_num  = diag(S'*K*S).^.5;   %Numerator
            se_denom = sum(Xtilde.^2);  %Denominator
            se_vec = (se_num./se_denom)*sqrt(T/(T-b_k));  %Standard errors
            pwr_stat_vec= (beta_hat_sim_cand(1,:)./se_vec')'; %tstats 
            %[t_val, se_kvar] = get_statistic(beta_hat_sim, eps_hat_sim, X_bs, h, l, s);
            %stat_vec(i) = t_val(1);
            %se_vec(i) = se_kvar(1);


%%  Get the empirical distribution of the statistic

            % pwr_stat_density = paretotails(pwr_stat_vec, 0.01, 0.99);
            pwr_stat_density = fitdist(abs(pwr_stat_vec), "Kernel", "Kernel", "normal");

            % [f, x1] = ecdf(pwr_stat_vec);
            % x = linspace(min(pwr_stat_vec), max(pwr_stat_vec), 100);
            % % type2_density = makedist('PiecewiseLinear', 'x', x1, 'Fx', f); 
            % ecdf(pwr_stat_vec);
            % hold on
            % plot(x,cdf(pwr_stat_density,x),'r--');
            % % hold on
            % % plot(x, cdf(type2_density,x),'b-');

%% Get the probability of rejection for the candidate beta
% using CV from testing the null given beta_cand =0, and the distribution
% of the statistic when testing H_0: beta ==0 given beta = beta_candidate
% get the power of the test
% Probability of rejection for the candidate beta

            % probs = cdf(pwr_stat_density, reshape(cv_array(jj,k,:),[1, 2])); % Probability of rejection for the candidate beta
            % prob_rejection = probs(1) + (1 - probs(2)); % Probability of rejection for the candidate beta
            % power_array(jj,k,kk) = prob_rejection;
            % fprintf('Probability of rejection for the candidate beta: %.4f\n', prob_rejection);

            probs = cdf(pwr_stat_density, cv_array(jj,k,1)); % Probability of rejection for the candidate beta
            prob_rejection = 1 - probs(1); % Probability of rejection for the candidate beta
            power_array(jj,k,kk) = prob_rejection;
        end 
    end
    pwrTimeEnd = toc(pwrTime);
    if verbose 
        fprintf('Pwr: Finished Simulations for %d Number of PCs in %.2f seconds \n ', sp, pwrTimeEnd); 
    end
    % fprintf('Finished Simulations for %d Number of PCs \n', sp);
% % % end
end

loopEnd = toc(loopStart);
if verbose 
    fprintf('Total time for simulations: %.2f seconds\n', loopEnd);
end
%% Collect results 

power_mean = mean(power_array,3); % Mean power for each l_cutoff and number of PCs

id_mat = zeros(length(l_cutoffs), length(PC_n));
[~, index] = max(power_mean,[],"all","linear");
id_mat(index) = 1; % Store the index of the maximum power

if size(power_mean,2) <= 1

    bw_tbl = array2table( [l_cutoffs' cv_array(:,:,1) cv_array(:,:,2) power_mean],'VariableNames', {'Cutoff','CV_Lower', 'CV_Upper', 'Mean Power'});
else
    bw_tbl = array2table( [l_cutoffs' power_mean],'VariableNames', ['Bwd/PC', string(PC_n) ]);
    cv_tbl = array2table( [l_cutoffs' cv_array(:,:,1)], 'VariableNames', ['Bwd/PC', string(PC_n) ]);
    % cv_u_tbl = array2table( [l_cutoffs' cv_array(:,:,2)], 'VariableNames', ['Bwd/PC', string(PC_n) ])
    % cv_l_tbl = array2table( [l_cutoffs' cv_array(:,:,1)], 'VariableNames', ['Bwd/PC', string(PC_n) ])
end

%% test using optimal CV, Bandwidth, and PCs

% CV
opt_cv = cv_array(index);
% bw 
bw_array = repmat(l_cutoffs',1,length(PC_n));
opt_bw = bw_array(index);

% PC
pc_array = repmat(PC_n,length(l_cutoffs),1);
opt_n_pc = pc_array(index);

% Run OLS
ii = min(opt_n_pc, size(W,2)); % Ensure sp does not exceed the number of PCs available
X_bs = [X W(:,1:ii)]; % X with B-splines, no intercept
[beta_hat, u_hat] = ols(y,X_bs,X_bs,'chol');
se_hat = kernel_var(u_hat,X_bs, X_bs, h, opt_bw, s, 'gaussian', 'chol', false);
t_stat = beta_hat(1)/se_hat(1);
rej_rule = abs(t_stat) > opt_cv;
ci_length = 2*abs(se_hat(1))*opt_cv;




%% Save the results

if save_results
    if verbose 
        fprintf('Saving results\n'); 
    end
    
    save(['outputs/' excrs 'bw_cv_sel.mat'])

    writetable(bw_tbl,['outputs/' excrs 'bw_cv_pwr.csv']);
    writetable(cv_tbl,['outputs/' excrs 'bw_cv.csv']);
    % writetable(cv_u_tbl,['outputs/' excrs 'bw_cv_u.csv']);
    % writetable(cv_l_tbl,['outputs/' excrs 'bw_cv_l.csv']);
end

%% Close parallel pool

if ~batchStartupOptionUsed
    if verbose 
        fprintf('Closing parallel pool\n'); 
    end
    delete(gcp('nocreate')) % Close the parallel pool
end


end