%% Setting directory and adding functions to path
cd Theta
addpath(genpath('./functions'))

%% Reading data

[y, X, s_x, s_y] =readvars('morgans_data.csv');

%% Running OLS and variance estimators

s=[s_x s_y];
Z=[ones(length(X),1) X];
[beta_hat, u_hat] = ols(y,Z,Z,'chol');
D_mat = getdistmat(s,0)
se_HR = HR_var(u_hat,Z,Z);
se_HAC = kernel_var(u_hat,Z,Z,D_mat,0.05,s,"uniform", "chol");
se_HAC_2 = kernel_var(u_hat,Z,Z,D_mat,0.1,s,"uniform", "chol");
se_HAC_3 = kernel_var(u_hat,Z,Z,D_mat,0.15,s,"uniform", "chol");
se_HAC_4 = kernel_var(u_hat,Z,Z,D_mat,0.2,s,"uniform", "chol");

%% Saving results

tbl=array2table([beta_hat se_HR se_HAC se_HAC_2 se_HAC_3 se_HAC_4],...
    'VariableNames', {'Coeffs' 'HR' 'Conley0.05' 'Conley0.10' 'Conley0.15' 'Conley0.20'})

writetable(tbl,'regs_tbl_morgans_data.csv');

%%  Using PC from Morgan

pc=readmatrix('morgans_splines_pc');
bam=readmatrix('morgans_splines_bam');
gam=readmatrix('morgans_splines_gam');

%% Running OLS
Z_pc=[ones(length(X),1) X pc];
[beta_hat_pc, u_hat_pc] = ols(y,Z_pc,Z_pc,'chol');

se_HR_pc = HR_var(u_hat_pc,Z_pc,Z_pc);
se_HAC_pc = kernel_var(u_hat_pc,Z_pc,Z_pc,D_mat,0.05,s,"uniform", "chol",true);
se_HAC_2_pc = kernel_var(u_hat_pc,Z_pc,Z_pc,D_mat,0.1,s,"uniform", "chol",true);
se_HAC_3_pc = kernel_var(u_hat_pc,Z_pc,Z_pc,D_mat,0.15,s,"uniform", "chol",true);
se_HAC_4_pc = kernel_var(u_hat_pc,Z_pc,Z_pc,D_mat,0.2,s,"uniform", "chol",true);


se_HAC_5_pc = kernel_var(u_hat_pc,Z_pc,Z_pc,D_mat,0.05,s,"uniform", "chol",false);
se_HAC_6_pc = kernel_var(u_hat_pc,Z_pc,Z_pc,D_mat,0.1,s,"uniform", "chol",false);
se_HAC_7_pc = kernel_var(u_hat_pc,Z_pc,Z_pc,D_mat,0.15,s,"uniform", "chol",false);
se_HAC_8_pc = kernel_var(u_hat_pc,Z_pc,Z_pc,D_mat,0.2,s,"uniform", "chol",false);

%% Saving results

tbl_pc=array2table([beta_hat_pc(1:2) se_HR_pc(1:2) se_HAC_pc(1:2) se_HAC_2_pc(1:2) se_HAC_3_pc(1:2) se_HAC_4_pc(1:2)],...
    'VariableNames', {'Coeffs' 'HR' 'Conley0.05' 'Conley0.10' 'Conley0.15' 'Conley0.20'})

writetable(tbl_pc,'regs_tbl_morgans_data_pc.csv');

tbl_pc_2=array2table([beta_hat_pc(1:2) se_HR_pc(1:2) se_HAC_5_pc(1:2) se_HAC_6_pc(1:2) se_HAC_7_pc(1:2) se_HAC_8_pc(1:2)],...
    'VariableNames', {'Coeffs' 'HR' 'Conley0.05' 'Conley0.10' 'Conley0.15' 'Conley0.20'})

writetable(tbl_pc_2,'regs_tbl_morgans_data_pc_not_fixed.csv');