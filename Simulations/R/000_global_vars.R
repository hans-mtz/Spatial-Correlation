# Setting global variables
# set.seed(66636)
n_obs = 250 #number of observations
n_sims = 1000 #number of simulations
rho_bar = 0.0301 # average pairwise correlation
beta = 0 # true beta of the model
fix_locations = TRUE
prod_dir = "Simulations/Products/"

save(list=ls(), file="Simulations/Products/glob_vars.RData")