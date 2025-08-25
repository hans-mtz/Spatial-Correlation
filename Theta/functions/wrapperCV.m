function [tstat_i, se_i] = wrapperCV(PC_n, W, X, eps_sim, beta_hat, h, s, b_cand, l_cutoffs, verbose)
% This function wraps the mainCVpar function to handle parallel processing
% and returns the t-statistics array for the given parameters.


% Initialize the output array
tstat_i = NaN(length(l_cutoffs), length(PC_n), length(b_cand));
se_i = NaN(length(l_cutoffs), length(PC_n), length(b_cand));
k = 0; % Initialize index for b_cand
    % Loop through each beta candidate, including the zero case
for b = b_cand
    loopStart = tic; % Start timer for critical values
    k = k + 1; % Increment index for beta candidates
    % j = 0; % Reset index for l_cutoffs
    j = 0; % Reset index for PC_n
    

    if verbose, fprintf('Starting Simulatios for beta_cand: %0.4f \n', b); end

    
    % Simulate data with the candidate beta and the simulated epsilon
    if (length(beta_hat) > 1)  
        beta_cand = [ b beta_hat(2:end)']; % Candidate beta and estimates of the Bsplines coefficients
    else
        beta_cand = b;
    end
    
    % [y_sim_pwr, eps_sim_pwr] = sim_w_tau(tau, h, X, beta_cand_pwr,B); % Drawing new eps_sim
    y_sim = X*beta_cand' + eps_sim; % Fixing eps_sim
    
    for sp = PC_n
        if verbose, fprintf('Starting Simulations for %d Number of PCs, and beta_cand = %.4f \n', sp, b); end
        pwrTime = tic; % Start timer for critical values
        ii = min(sp, size(W,2)); % Ensure sp does not exceed the number of PCs available
        X_bs = [X W(:,1:ii)]; % X with B-splines, no intercept
        [beta_hat_sim_cand, eps_hat_sim_cand] = ols(y_sim, X_bs, X_bs, 'chol');
        
        j = j + 1; % Increment index for l_cutoffs
        i = 0; % Reset index for l_cutoffs
        for l = l_cutoffs %variance changes per l_cutoff
            i = i + 1; % Increment index for l_cutoffs
            % For CV Test H_0: beta_0 = 0 | beta = beta_cand
            [t_val, se_kvar] = get_statistic(beta_hat_sim_cand, eps_hat_sim_cand, X_bs, h, l, s);
            tstat_i(i,j,k) = t_val(1);
            se_i(i,j,k) = se_kvar(1);
        end 
        pwrTimeEnd = toc(pwrTime);
        if verbose, fprintf('Finished Simulations for %d Number of PCs in %.2f seconds \n ', sp, pwrTimeEnd); end
    end
    if verbose, fprintf('Finished Simulations for beta_cand: %0.4f in %0.2f seconds \n', b, toc(loopStart)); end
end
 
end