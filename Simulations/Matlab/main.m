% %-----------------------------
% % Objective: This code generates a random sample of spatially correlated
% % data and computes the kernel estimator of the OLS variance a la
% % Conely(1999).
% % Author: Hans Martinez hmarti33@uwo.ca
% % Date: July 21, 2023
% % ----------------------------
clear;
if batchStartupOptionUsed
  addpath(genpath('./SCPC_Matlab_Implementation_20211123'))
  addpath(genpath('./21200057'))
end

% Table 1

excercise='light_locs_me_3p_matching_rhos_tims_stat';
rho_bar_true=0.03; %DGP's rho; MW's rho=0.03
% spatial = 0;% 1: Uniform dist locs: other: Hihgly concentrated locs
% e=0.03; %Measurement error
run('simulate_optimal_me.m')
clear;

% Table 2
% excercise='light_locs_me_3p_true_rho_above_tims_stat';
% rho_bar_true=0.06; %DGP's rho; MW's rho=0.03
% run('simulate_optimal_me.m')
% clear;
% % Table 3
% excercise='light_locs_me_3p_true_rho_below_tims_stat';
% rho_bar_true=0.01; %DGP's rho; MW's rho=0.03
% run('simulate_optimal_me.m')
% clear;





