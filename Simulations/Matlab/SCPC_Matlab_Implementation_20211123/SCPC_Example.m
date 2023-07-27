% This is main driver program for SCPC .. as you will see the program
% (a) Reads in the dataset
% (b) Runs OLS regression and computes White standard errors
% (c) Runs SCPC to produce SC-robust SEs and critical values
%
% A 95% SC-robust confidence interval is provided
% Also, critical values are listed for 68, 90, 95 and 99 percent confidence intervals, which can be formed as the OLS estimate plus and minus the SC-robust SE times the critical value
% 
% Two sample data sets are provided. The first has small n (n= 74) and the second has large n (n = 10,000).
% 
% Note that the SCPC analysis uses rhobar ... a user-supplied upper bound on the average spatial correlation of the relevant observations. The default value of 0.03 is used here. 

clear all;                
this_date = datestr(now,'yyyymmdd');

%%%%%%%%%%%%%%%%%%%%% (a) Read in the data set %%%%%%%%%%%%%%%%%%

% Read in the Data used for this test
% Read a small data set  -- These are from STATA test data set
T = readtable('scpc_testdata_auto.xlsx');
mpg = table2array(T(:,1));
weight = table2array(T(:,2));
length = table2array(T(:,3));
s = table2array(T(:,end-1:end));
n = size(s,1);

%{
% Read a large data -- these are random numbers, but with same form as scpc_testdata_auto.xlsx
% Note .. this takes approximately 15 seconds to execute on my computer .. it may take longer
% on yours
T = readtable('large_data.xlsx');
mpg = table2array(T(:,1));
weight = table2array(T(:,2));
length = table2array(T(:,3));
s = table2array(T(:,end-1:end));
n = size(s,1);
%}

tic;

%%%%%%%%%%%%%%%%%%%%% (b) OLS and White Standard Errors %%%%%%%%%%%%%%%%%%
% Set up Regression
Y = mpg;
X = [weight length ones(n,1)];
k = size(X,2);

% Form XX, XY, etc;
XX = (X'*X);
XY = (X'*Y);
XX_inv = inv(XX);
S_XX = (X'*X)/n;
S_XX_inv = inv(S_XX);

% OLS
beta_hat = XX_inv*XY;

% Residuals
u_hat = Y-X*beta_hat;
Xu = X.*repmat(u_hat,1,k);  % X*u

% White Covariance
Xu_Xu = Xu'*Xu;
V_beta = XX_inv*Xu_Xu*XX_inv';
% STATA DF adjustment
V_beta = V_beta*n/(n-k);
SE_beta_hat = sqrt(diag(V_beta));

%%%%%%%%%%%%%%%%%%%%% (c) SCPC results %%%%%%%%%%%%%%%%%%
y_data = repmat(beta_hat',n,1)+Xu*S_XX_inv';     % j'th column is beta_hat_j+appropiate linear combination of Xu
rhobar = 0.03;     % Upper bound on average spatial correlation
ci_level = 0.95;   % Confidence level
latlongflag = 0;   % 0 if Euclidean distance; 1 if s = [lattitude longitude]
nobs_s = size(s,1);
s1 = unique(s,'rows');
nobs_unique = size(s1,1);

rslt = scpc(y_data,s,latlongflag,rhobar,ci_level);  % Run SCPC

% Print Out Some Results
% OLS and White SE Robusts
fprintf('OLS Results with iid (Heter-robust, White) SEs a la STATA \n');
fprintf('Regressor      Coef          SE  \n');
fprintf('weight       %5.4G   %5.4g  \n',[beta_hat(1) SE_beta_hat(1)]);
fprintf('length       %5.4G   %5.4g  \n',[beta_hat(2) SE_beta_hat(2)]);
fprintf('_const       %5.4G   %5.4g  \n',[beta_hat(3) SE_beta_hat(3)]);

fprintf('\n\n');
fprintf('Results from SCPC \n');
if latlongflag == 0
    fprintf('Euclidean norm used to compute distance between locations is used \n')
elseif latlongflag == 1
    fprintf('Locations are lattitude and longitude: using distances on a sphere \n')
end
fprintf(['Number of observations:' num2str(nobs_s) '\n']);
fprintf(['Number of unique locations:' num2str(nobs_unique) '\n']);
fprintf('Spatial correlation robust SEs and confidence intervals for rhobar = %5.4f \n',rhobar);
fprintf('Optimal value of q = %3i for minimum length 95 percent conf. interval \n',rslt.q);

% OLS and SCPC
fprintf('OLS Results with SCPC robust SEs and 95 percent CI \n');
fprintf('Regressor      Coef          SE    Conf Int() \n');
fprintf('weight       %5.4G   %5.4g  (%5.4g,%5.4g)  \n',[rslt.beta_hat(1) rslt.se_beta_hat(1) rslt.ci(1,:)]);
fprintf('length       %5.4G   %5.4g  (%5.4g,%5.4g)\n',[rslt.beta_hat(2) rslt.se_beta_hat(2) rslt.ci(2,:)]);
fprintf('_const       %5.4G   %5.4g  (%5.4g,%5.4g)\n',[rslt.beta_hat(3) rslt.se_beta_hat(3) rslt.ci(3,:)]);
fprintf('\n\n');
fprintf('Critical values for tests or CIs');
fprintf('Level (CI Coverage) Critical Value: \n');
for i = 1:size(rslt.cv_vec,1);
    fprintf('%5.3f  (%5.3f):  %5.4f \n',[1-rslt.ci_level_vec(i) rslt.ci_level_vec(i) rslt.cv_vec(i)]);
end;

toc;
    