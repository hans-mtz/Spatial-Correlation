%-----------------------------
% Objective: This code generates a random sample of spatially correlated
% data and computes the kernel estimator of the OLS variance a la
% Conely(1999).
% Author: Hans Martinez hmarti33@uwo.ca
% Date: July 21, 2023
% ----------------------------
%----- Generating data ------%
% Setting up parameters
global T k beta rho_bar b L
T = 250;
k= 3;
beta = [ 1; 1.5; 0.5];
rho_bar = 0.03;
b = 0.08;
L = 0.08;


% Generate variables
s = rand(T,1); % vector of locations
X = [ones(T,1) randn(T,k-1)]; % matrix of independent variables (can be endogenous)
Z = randn(T,k); % matrix of exogenous variables
D_mat = getdistmat(s,false); % matrix of distances
c = getcbar(rho_bar, D_mat); % calculating the c_min given rho_bar
Sigma = exp(-c*D_mat); % var-cov matrix 
rho_estimated = mean(reshape(Sigma,[],1)); % estimated c_min corresponds to 
fprintf('Estimated c_min corresponds to a rho_bar estimated of %5.4g \n', rho_estimated);
mu = zeros(T,1); % mean zero for error terms
u = mvnrnd(mu,Sigma)'; % sampling from G_exp(c_min)
y = X*beta + u;

% Saving data to open in Stata -----------

data = table(y, X, s);
writetable(data, 'Stata/data.csv');

%---------OLS-------------%
XX = X'*X;
ZX = Z'*X;
Zy = Z'*y;
Xy = X'*y;
ZX_inv = inv(ZX);
XX_inv = inv(XX);

beta_IV_hat = ZX_inv*Zy;
beta_hat = XX_inv*Xy;
u_hat = y - X*beta_hat;
fprintf('Checking moment condition E[Zu]=0 %5.4g \n', Z'*u_hat/T);
fprintf('Checking moment condition E[Xu]=0 %5.4g \n', X'*u_hat/T);


%------Variance kernel estimator------------%
% Zu_hat = Z.*repmat(u_hat,1,k);
Xu_hat = X.*repmat(u_hat,1,k);

V_hat = zeros(k); % Initializing V matrix

for i=1:T
%     for j=find(D_mat(i,:)<=L)
      for j=i
        V_hat = V_hat + Xu_hat(i,:)'*Xu_hat(j,:);
    end
end

% V_hat0 = V_hat/T; No need to make comparable to Stata's HR SE 
% when K=1 if i=j; 0 otherwise

% V_kernel = (XX_inv)*V_hat*(XX_inv)'*T/(T-k); Stata's DF
V_kernel = (XX_inv)*V_hat*(XX_inv)'/T;
SE = sqrt(diag(V_kernel));

%-------Stata's HR ------------------------%
V_HR0 = Xu_hat'*Xu_hat;
V_HR = XX_inv*V_HR0*XX_inv';
V_HR = V_HR*T/(T-k);
SE_HR = sqrt(diag(V_HR));






