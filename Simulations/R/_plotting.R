load("Simulations/Products/functions.RData")

# Plot covariance vs distance ---------

s_l = runif(n_obs) #vector of locations
d_mat <- abs(outer(s_l,s_l,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
c_min <- get_c(d_mat,rho_bar)
Sigma_c <- exp(-c_min*d_mat)
df <- model_1()

df_plot <- stack(data.frame(d_mat)) |> cbind(stack(data.frame(Sigma_c)))
names(df_plot) <- c("distance","ind","cov","ind2")

df_plot$s_1 <- s_l

with(
    df_plot,
    plot(distance,cov, main="Distancs vs. Covariance")
)

st(df_plot[,c("distance","cov","s_1")])


# Plot hist of covariances ---------
df_plot
five_num<-with(
    df_plot,
    fivenum(cov)
    # boxplot.stats(cov)
)
names(five_num)<-c("Min","1st Q","Median","3rd Q.","Max")

mean(df_plot$cov, na.rm=FALSE)

with(
    df_plot[df_plot$distance!=0 & df_plot$cov<1 & df_plot$cov>0.01,],
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

save(df_plot, five_num, file=paste0(prod_dir,"plots.RData"))
