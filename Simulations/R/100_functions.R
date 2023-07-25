# library(RStata)
# options("RStata.StataPath"="/Applications/Stata/StataSE.app/Contents/MacOS/stata-se")
# options("RStata.StataVersion" = 16)
load("Simulations/Products/glob_vars.RData")


## Calibrating the critical value ----

mean_pairwise_corr <- function(Sigma){
    sum(Sigma)/(nrow(Sigma)*(ncol(Sigma)-1))
}

rho_c <- function(c,d){
    sigma <- exp(-c*d)
    rho <- mean_pairwise_corr(sigma)
    return(rho)
}

obj_f <- function(x,r,d) rho_c(x,d)-r

get_c <- function(d,r){
  ans<-uniroot(obj_f,c(0.1,1000),tol=1e-8,r=r,d=d)
  return(ans$root[1])
}

if (fix_locations) {
    print("Fixing locations")
    s_l = runif(n_obs) #vector of locations
    d_mat <- abs(outer(s_l,s_l,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
    c_min <- get_c(d_mat,rho_bar)
}

model_1 <- function(n=n_obs, c=c_min,d=d_mat, s=s_l) {
    # n <- 250 #n of observations
    # n_sim <- 100 #n of simulations
    #c <- 75.167 #critical value of the covariance function
    # s_1 <- runif(n) #location vector
    # d <- abs(outer(s_1,s_1,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
    # c <- get_c(d,r)
    Sigma_c <- exp(-c*d) # Covariance function -> COV matrix
    e_l <- mvtnorm::rmvnorm(1,sigma=Sigma_c) # sampling from gaussian process
    y <- beta*1+ t(e_l)
    df <- data.frame(y=y,x1=1,e_l=t(e_l), s_1=s)
    return(df)
}

model_2 <- function(n=n_obs, c=c_min,d=d_mat, s=s_l) {
    # n <- 250 #n of observations
    # n_sim <- 100 #n of simulations
    #c <- 75.167 #critical value of the covariance function
    # s_1 <- runif(n) #location vector
    # d <- abs(outer(s_1,s_1,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
    # c <- get_c(d,r)
    c_2 <- c/2
    Sigma_c_2 <- exp(-c_2*d) # Covariance function -> COV matrix
    e_l <- mvtnorm::rmvnorm(1,sigma=Sigma_c_2) # sampling from gaussian process
    x_l <- mvtnorm::rmvnorm(1,sigma=Sigma_c_2)
    # if (!length(e_l)==length(x_l)) print("e_l and x_l are not ok")
    y <- beta*t(x_l)+ t(e_l)
    df <- data.frame(y=y,x1=t(x_l),e_l=t(e_l), s_1=s)
    return(df)
}

model_3 <- function(n=n_obs, c=c_min,d=d_mat, s=s_l) {
    # n <- 250 #n of observations
    # n_sim <- 100 #n of simulations
    #c <- 75.167 #critical value of the covariance function
    # c <- c/2
    # s_1 <- runif(n) #location vector
    # d <- abs(outer(s_1,s_1,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
    # c <- get_c(d,r)
    Sigma_c <- exp(-c*d) # Covariance function -> COV matrix
    e_l <- mvtnorm::rmvnorm(1,sigma=Sigma_c) # sampling from gaussian process
    cutoff <- quantile(s, probs = 0.85)[1]
    x_l <- ifelse(s>cutoff,0.85,-0.15) #Changed to <= instead of <
    y <- beta*x_l+ t(e_l)
    df <- data.frame(y=y,x1=x_l,e_l=t(e_l), s_1=s)
    return(df)
}

model_4 <- function(n=n_obs, c=c_min,d=d_mat, s=s_l) {
    # n <- 250 #n of observations
    # n_sim <- 100 #n of simulations
    #c <- 75.167 #critical value of the covariance function
    # c <- c/2
    # s_1 <- runif(n) #location vector
    # d <- abs(outer(s_1,s_1,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
    # c <- get_c(d,r) # Get c value 
    Sigma_c <- exp(-c*d) # Covariance function -> COV matrix
    e_l <- mvtnorm::rmvnorm(1,sigma=Sigma_c) # sampling from gaussian process
    # cutoff <- quantile(s_l, probs = 0.85)[1]
    x_l <- cumsum(c(0,rnorm(n,0,rw_sigma))) # Random walk with mean zero
    # x_l_mean <- mean(x_l)
    # demeaned_rw <- x_l[2:(n+1)] - x_l_mean
    # y <- beta*demeaned_rw + t(e_l)

    y <- beta*x_l[2:(n+1)]+ t(e_l) # Changed from demeaned RW
    df <- data.frame(y=y,x1=x_l[2:(n+1)],e_l=t(e_l), s_1=s)
    return(df)
}

ratio_list<-function(x,f=model_1){
    if (!fix_locations) {
        print("Randomizing locations")
        s_l <- runif(n_obs) #vector of locations
        d_mat <- abs(outer(s_l,s_l,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
        c_min <- get_c(d_mat,rho_bar)
        df <- f(c=c_min,d=d_mat,s=s_l)
    } else {
        df <- f()
    }
    
    stata_df <- stata(
        "Simulations/Stata/cmd_prog.do", 
        data.in=df, data.out=TRUE
    )
    ## Using p-values
    # p_val_HR <-min(stata_df$p_val_HR)[1]
    # p_val_SCPC<-min(stata_df$p_val_SCPC)[1]
    # ratio_HR <- ifelse(p_val_HR>0.05,1,0)
    # ratio_SCPC <- ifelse(p_val_SCPC>0.05,1,0)

    ## Using t-statistic
    t_HR <-min(stata_df$t)[1]
    t_SCPC<-min(stata_df$t_SCPC)[1]
    t_SCPC_cv<-min(stata_df$t_SCPC_cv)[1]

    ratio_HR <- ifelse(abs(t_HR)>1.960,1,0)
    ratio_SCPC <- ifelse(abs(t_SCPC)>abs(t_SCPC_cv),1,0)

    # if(pval>0.05){
    #     ratio = 1
    # }else{
    #     ratio = 0
    # }
    return(data.frame(HR=ratio_HR,SCPC=ratio_SCPC))
}

rejection_ratio<-function(f=model_1,nsims=n_sims){
    print("Started working on another model")
    # if (nsims<=500) {
    #     r_list<-lapply(1:nsims, ratio_list,f=f)
    #     r_df<-do.call(rbind,r_list)
    #     r <- colSums(r_df)/nsims
    # } else {
    #    chunk <- c(seq(1,nsims,500),nsims)
    #    r_chunks <- list()
    #    for (i in 1:(length(chunk)-1)) {
    #         r_list<-lapply(chunk[i]:(chunk[i+1]-1), ratio_list,f=f)
    #         r_df<-do.call(rbind,r_list)
    #         r_chunks[i] <- colSums(r_df)/(chunk[i+1]-chunk[i])
    #         save(
    #             r_chunks,
    #             file= paste(
    #                 "Simulations/Products/ratio",
    #                 as.character(quote(f())),
    #                 "n_chunk",
    #                 (chunk[i+1]-1),
    #                 ".RData",
    #                 sep = "_"
    #             )
    #         )
    #     }
    #     r <- colSums(do.call(rbind,r_chunks))/length(r_chunks)
    # }
    r_list<-lapply(1:nsims, ratio_list,f=f)
    r_df<-do.call(rbind,r_list)
    r <- colSums(r_df)/nsims
    # print(as.character(quote(f())))
    save(
        r,
        file= paste(
            "Simulations/Products/ratio",
            as.character(quote(f())),
            "nsim",
            n_sims,
            ".RData",
            sep = "_"
        )
    )
    return(r)
}

save(list=ls(), file=paste0(prod_dir,"functions.RData"))