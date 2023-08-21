
ts_df <- read.csv("Simulations/Products/data_ts.csv")
u_hat <- residuals(lm(y~X_2,data=ts_df))
lm(y~X_2,data=ts_df) |> summary()
lm(y~X_2,data=ts_df) |> BIC()

my_BIC <- function(resid,k,method){
    n<-length(resid)
    SSR<- sum(resid^2)
    sigma_2 <- var(resid)
    if (method=="SSR") {
        return(n*log(SSR)+(k+1)*log(n)-n*log(n))
    } else {
       return(n*log(sigma_2)+(k+1)*log(n))
    }
}


my_BIC(u_hat,2,"SSR")
my_BIC(u_hat,2,"")
lm(y~X_2,data=ts_df) |> BIC()
lm(y~X_2,data=ts_df) |> AIC()