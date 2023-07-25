# Setting global variables
# set.seed(66636)
n_obs = 250 #number of observations
n_sims = 10000 #number of simulations
rho_bar = 0.0301 # average pairwise correlation
beta = 0 # true beta of the model
fix_locations = TRUE
prod_dir = "Simulations/Products/"
rw_sigma = 3 #Variance of the random walk process

save(list=ls(), file=paste0(prod_dir,"glob_vars.RData"))