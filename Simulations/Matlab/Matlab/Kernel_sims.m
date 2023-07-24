% Simulations - Rejection Frequencies for Kernel Var Estimator%
N_models = 2;
N_reps = 200;

% Setting up parameters
global T k beta rho_bar b L
T = 250;
k= 3;
beta = [ 1; 1.5; 0.5];
rho_bar = 0.03;
b = 0.08;
L = 0.001;

% Fixing locations
s = rand(T,1); % vector of locations
rho = [rho_bar, 0.0];
rej_freq = zeros(N_models,k);
for m=1:N_models
%     rej_freq_m = 0;
    rej_freq_m = zeros(1,k);

    for r=1:N_reps
        % generate data      
        [y, X, D_mat] = DGP(beta,T,k,s,rho(m),1);
%         writetable(table(y, X, s), strcat('../Stata/data_',num2str(m),'.csv')); % Saving
        [beta_hat, SE] = kernel_var(y,X,D_mat,k,L,T);
        t_stat = (beta_hat - beta)./SE;
        rej_i = abs(t_stat) > 1.96;
        rej_freq_m = rej_freq_m + rej_i';      
    end
    rej_freq(m,:) = rej_freq_m./N_reps;
end

fprintf('Distance cutoff; %5.4g\n', L);
fprintf('                                 cons  beta_1 beta_2\n')
fprintf('Rejection frequency spatial was (%5.4g %5.4g %5.4g)\n',rej_freq(1,:));
fprintf('Rejection frequency iid was     (%5.4g %5.4g %5.4g)\n',rej_freq(2,:));
        
        