function [sims_stats_mat,raw_data_bs,raw_data_bs_ci,n_imaginary] = wrapper_PC(rho,l_cutoff,beta,theta,s, W, fix,h)

    raw_data_bs =NaN(size(W,2),length(beta));
    n_imaginary =zeros(size(W,2),length(beta));
    sims_stats_mat =NaN(size(W,2),4);

    % CI
    raw_data_bs_ci =NaN(size(W,2),length(beta));

    % Generate data 

    [y, X, ~] = DGP(theta,s,rho,false,h);

    %% Kernel HAC with PC of B-Splines 

    for j=1:size(W,2)
        
        S = W(:,1:j);

        X_bs =[X S]; % Removing intercept in GDP

        %% OLS with B-Splines

        [beta_hat_bs, u_hat_bs] =ols(y,X_bs,X_bs,'chol');
        sims_stats_mat(j,3) =get_nn_corr(u_hat_bs,h,1);
        sims_stats_mat(j,2) =bic(u_hat_bs,X_bs,'hansen');

        % Null Hypothesis for OLS with B-Splines

        H_0_bs = beta_hat_bs(1:length(beta)) - beta';

        % Testing with HR

        se_hr_bs = HR_var(u_hat_bs,X_bs,X_bs);

        sims_stats_mat(j,4) = se_hr_bs(length(beta));
        t_stat_hr_bs = H_0_bs./se_hr_bs(1:length(beta));
        rej_rule_hr_bs = abs(t_stat_hr_bs) > 1.96;
        sims_stats_mat(j,1) = rej_rule_hr_bs(length(beta));


            
        %% t-Tests %%
        %%%% Kernel HAC %%%%

        %% Triangle Kernel HAC
        se_bs = kernel_var(u_hat_bs,X_bs,X_bs,h,l_cutoff,s,'gaussian','chol',fix);

        % Triangle Kernel

        %%%% Kernel HAC + B-SPlines %%%%
        % t statistic and rejection for kernel+bsplines  

        % Triangle Kernels
        n_imaginary(j,:) = any(se_bs~=real(se_bs));
        for k=1:length(beta)
            if (se_bs(k)==real(se_bs(k)))
                t_stat_bs = H_0_bs(k)/se_bs(k);
                rej_rule_bs = abs(t_stat_bs) > 1.96;
                raw_data_bs(j,k) = rej_rule_bs;
                raw_data_bs_ci(j,k) = 2.*abs(se_bs(k)).*1.96; 
            end
        end

    end

end