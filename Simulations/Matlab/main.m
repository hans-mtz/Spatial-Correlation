% %-----------------------------
% % Objective: This code generates a random sample of spatially correlated
% % data and computes the kernel estimator of the OLS variance a la
% % Conely(1999).
% % Author: Hans Martinez hmarti33@uwo.ca
% % Date: July 21, 2023
% % ----------------------------
if batchStartupOptionUsed
  addpath(genpath('./SCPC_Matlab_Implementation_20211123'))
  % addpath(genpath('./21200057'))
end

% run('crank_bs_sims_chol')
% run('crank_bs_ear1_sims')
% run('grid_sims')
% run('e_ar1')
run('bic_spline_sims')
run('plot_nspline_hist_by_min_bic')







