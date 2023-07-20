load("Simulations/Products/functions.RData")

# Plot covariance vs distance
s_l = runif(n_obs) #vector of locations
d_mat <- abs(outer(s_l,s_l,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
c_min <- get_c(d_mat,rho_bar)
Sigma_c <- exp(-c_min*d_mat)
df <- model_1()

df_plot <- stack(data.frame(d_mat)) |> cbind(stack(data.frame(Sigma_c)))
names(df_plot) <- c("distance","ind","cov","ind2")

dist_cov_plot <- with(
    df_plot,
    plot(distance,cov, main="Distancs vs. Covariance")
)

# Plot hist of covariances

with(
    df_plot,
    # hist(distance,breaks=100, freq=TRUE)
    hist(cov, breaks=50,freq=T, main="Covariance histogram")
)

with(
    df_plot,
    hist(distance,freq=TRUE, main="Distance histogram")
    # hist(cov, breaks=300, freq=T)
)

with(
    df,
    hist(s_1,freq=TRUE, main="Locations histogram")
    # hist(cov, breaks=300, freq=T)
)

save(list=grep("plot",ls(), value=TRUE), file=paste0(prod_dir,"plots.RData"))
