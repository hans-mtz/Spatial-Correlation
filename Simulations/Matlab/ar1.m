%% AR1 residualts %%

% Setting up parameters
global T k beta rho_bar L
T = 250; % obs
beta = [ 1; 1.5; 0.5];
k= length(beta); % number of covariates
rho_bar = 0.03; % Average pairwise correlation
L = [0.03 0.05]; % Vector of cutoff distances
R = [64 32 12 5]; % Number of splines
cholesky_flag = 'ldl';
% N_order_splines = 2;
% N_models = 2;
% N_reps = 250;
% 
% N_cutoffs = length(L);
% N_triangles = length(R);
% N_estimators = 2+N_cutoffs+N_cutoffs*N_triangles*N_order_splines; 


% Fixing locations
s = rand(T,1); % vector of locations

% Generating data

[y, X, D_mat] = DGP(beta, k, s, rho_bar, 2);

% Getting B-SPline matrix

S=get_bspline_mat(s,4,2);

[s_srtd, sorting_index] = sort(s); % Sorting by locations

% OLS and residuals

[beta_hat, u_hat] = ols(y, [X(:, 2:end) S], [X(:, 2:end) S], cholesky_flag);

%% AR1 on residuals

u_hat_sorted = u_hat(sorting_index);

u_hat_t = u_hat_sorted(1:end-1);

u_hat_t_1 =  u_hat_sorted(2:end);

[beta_u_hat] = ols(u_hat_t, [ones(T-1,1) u_hat_t_1], [ones(T-1,1) u_hat_t_1], cholesky_flag);
    

%% %%

fprintf('ModeL           Coef. beta1  beta2 B-Splines slopes \n');
fprintf('Kernel Splines        %5.4g %5.4g  %5.4g %5.4g %5.4g %5.4g \n', beta_hat);
fprintf('AR1 residuals   %5.4g %5.4g \n', beta_u_hat); 

%%

Sigma = diag(ones(T,1))+diag(ones(T-1,1).*0.79,-1);

try chol(Sigma)
    disp('symmetric positive definite')
catch ME
    disp('Not PSD')
end
