function [sims_stats_mat,raw_data_k_u,raw_data_k_u_ci,raw_data_k,raw_data_k_ci,raw_data_bs,raw_data_bs_ci,raw_data_bs_u,raw_data_bs_u_ci,n_imaginary,n_imaginary_u] = wrapper_opt_PC(rho,l_cutoffs,beta,theta,s, W, n_dropped, fix,h)

    raw_data_k_u=NaN(length(rho),length(l_cutoffs),length(beta));
    raw_data_bs_u=NaN(length(rho),length(l_cutoffs),length(beta));
    n_imaginary_u=zeros(length(rho),length(l_cutoffs),length(beta));

    raw_data_k=NaN(length(rho),length(l_cutoffs),length(beta));
    raw_data_bs=NaN(length(rho),length(l_cutoffs),length(beta));
    n_imaginary=zeros(length(rho),length(l_cutoffs),length(beta));

    sims_stats_mat=NaN(length(rho),6,2);

    % CI
    raw_data_k_u_ci=NaN(length(rho),length(l_cutoffs),length(beta));
    raw_data_bs_u_ci=NaN(length(rho),length(l_cutoffs),length(beta));

    raw_data_k_ci=NaN(length(rho),length(l_cutoffs),length(beta));
    raw_data_bs_ci=NaN(length(rho),length(l_cutoffs),length(beta));

    sims_stats_mat(:,2,2) = n_dropped;

    for j=1:length(rho)

        % Generate data for a rho
        [y, X, ~] = DGP(theta,s,rho(j),false,h);

        %% OLS - 
        [beta_hat, u_hat] = ols(y,X,X,'chol');
        sims_stats_mat(j,1,1) = get_nn_corr(u_hat,h,1);
        sims_stats_mat(j,4,1) = bic(u_hat,X,'hansen');

        % Null hypothesis
        H_0 = beta_hat - beta';

        % Testing with HR
        se_hr = HR_var(u_hat,X,X);
        

        % t_stat = H_0./se_hr;
        % rej_rule_hr = abs(t_stat) > 1.96;
        % sims_stats_mat(j,3,1) = rej_rule_hr;

        % Ignoring imaginary se

        sims_stats_mat(j,6,1) = se_hr(length(beta));
        t_stat_hr = H_0/se_hr;
        rej_rule_hr = abs(t_stat_hr) > 1.96;
        sims_stats_mat(j,3,1) = rej_rule_hr(length(beta));



        %% Kernel HAC with B-Splines %%
        
        % [n_splines, fval] = get_opt_splines_grid(s,splines_order,y,X,'nn','chol'); % Find optimal number of splines, Step functions
        % [S, ~, n_dropped] = get_bsplines(s,n_splines,splines_order); %Fixed number of splines
        % S = readmatrix('morgans_splines_pc'); n_dropped =0; %Use Morgans's Splines
        % sims_stats_mat(j,2,2) = n_dropped;
        % sims_stats_mat(j,5,2) = n_splines;
        % X_bs
        % disp(size(X))
        % disp(size(S))

        %% Get optimal number of PCs of B-Splines

        [n_PCs, ~] = get_opt_n_pcs_grid(y,X,W,h,'nn','chol');
        % disp(n_PCs)
        X_bs = [X W(:, 1:n_PCs)]; % Removing intercept in GDP
        sims_stats_mat(j,5,2) = n_PCs;
        %% OLS with B-Splines
        [beta_hat_bs, u_hat_bs] = ols(y,X_bs,X_bs,'chol');
        sims_stats_mat(j,1,2) = get_nn_corr(u_hat_bs,h,1);
        sims_stats_mat(j,4,2) = bic(u_hat_bs,X_bs,'hansen');

        % Null Hypothesis for OLS with B-Splines
        H_0_bs = beta_hat_bs(1:length(beta)) - beta';

        % Testing with HR
        se_hr_bs = HR_var(u_hat_bs,X_bs,X_bs);
        % t_stat_hr_bs = H_0./se_hr_bs(1:length(beta));
        % rej_rule_hr_bs = abs(t_stat_hr_bs) > 1.96;
        % sims_stats_mat(j,3,2) = rej_rule_hr_bs;

        % Ignoring imaginary se

        sims_stats_mat(j,6,2) = se_hr_bs(length(beta));
        t_stat_hr_bs = H_0_bs./se_hr_bs(1:length(beta));
        rej_rule_hr_bs = abs(t_stat_hr_bs) > 1.96;
        sims_stats_mat(j,3,2) = rej_rule_hr_bs(length(beta));


        % i=1;
        % for l=l_cutoffs
        for i=1:length(l_cutoffs)
            
            %% t-Tests %%
            %%%% Kernel HAC %%%%
            % SEs - varies by cutoff
            %% Uniform Kernel HAC
            se_u = kernel_var(u_hat,X,X,h,l_cutoffs(i),s,'uniform','chol',fix);
            se_bs_u = kernel_var(u_hat_bs,X_bs,X_bs,h,l_cutoffs(i),s,'uniform','chol',fix);
            %% Triangle Kernel HAC
            se = kernel_var(u_hat,X,X,h,l_cutoffs(i),s,'gaussian','chol',fix);
            se_bs = kernel_var(u_hat_bs,X_bs,X_bs,h,l_cutoffs(i),s,'gaussian','chol',fix);

            % Uniform Kernel
            % if (se_u==real(se_u))
            t_stat_u = H_0./se_u;
            rej_rule_u = abs(t_stat_u) > 1.96;
            raw_data_k_u(j,i,:) = rej_rule_u;
            raw_data_k_u_ci(j,i,:) = 2.*abs(se_u).*1.96;
            % end

            % Triangle Kernel
            % if (se==real(se) )
            t_stat = H_0./se;
            rej_rule = abs(t_stat) > 1.96;
            raw_data_k(j,i,:) = rej_rule;
            raw_data_k_ci(j,i,:) =2.*abs(se).*1.96;
            % end

            %%%% Kernel HAC + B-SPlines %%%%
            % t statistic and rejection for kernel+bsplines  

            % Triangle Kernels
            n_imaginary(j,i,:) = any(se_bs~=real(se_bs));
            for k=1:length(beta)
                if (se_bs(k)==real(se_bs(k)))
                    t_stat_bs = H_0_bs(k)/se_bs(k);
                    rej_rule_bs = abs(t_stat_bs) > 1.96;
                    raw_data_bs(j,i,k) = rej_rule_bs;
                    raw_data_bs_ci(j,i,k) = 2.*abs(se_bs(k)).*1.96; 
                    % error('se are imaginary in Kernel + Bsplines');
                else
                    % raw_data_bs(j,i,k) = NaN;
                    % n_imaginary(j,i,k)=1;
                end
            end

            % Uniform Kernels
            n_imaginary_u(j,i,:)=any(se_bs_u~=real(se_bs_u));
            for k=1:length(beta)
                if (se_bs_u(k)==real(se_bs_u(k)))
                    t_stat_bs_u = H_0_bs(k)/se_bs_u(k);
                    rej_rule_bs_u = abs(t_stat_bs_u) > 1.96;
                    raw_data_bs_u(j,i,k) = rej_rule_bs_u;
                    raw_data_bs_u_ci(j,i,k) = 2.*abs(se_bs_u(k)).*1.96; 
                    % error('se are imaginary in Kernel + Bsplines');
                else
                    % raw_data_bs_u(j,i,k) = NaN;
                    % n_imaginary_u(j,i,k)=1;
                end
            end

            % i=i+1;
        end
        % j=j+1;
    end
end