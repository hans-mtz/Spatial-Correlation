# GDP for simulations -----
n <- 250 #n of observations
n_sim <- 100 #n of simulations
c <- 75.167 #critical value of the covariance function
s_1 <- runif(n) #location vector
# s_l <- lapply(seq(1,n_sim,by=1),(\(x)runif(n)))
# S_l <- do.call(cbind, s_l)
d <- abs(outer(s_1,s_1,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
Sigma_c <- exp(-c*d) # Covariance function -> COV matrix
e_l <- mvtnorm::rmvnorm(1,sigma=Sigma_c) # sampling from gaussian process
y <- 1+ t(e_l)
df <- data.frame(y=y,x1=1,e_l=t(e_l), s_1=s_1)

stata("")

## Calibrating the critical value ----

lapply(seq(75.16,75.17,by=0.001),(\(x)rho_c(x,d)))
mean_pairwise_corr(Sigma_c)

## Generating DF -----

example <- data.frame(y=y,x1=1,e_l=t(e_l), s_1=s_1)
example <- model_1()
write.csv(example, file="Simulations/example.csv")

## Stata in R -----

library(RStata)
options("RStata.StataPath"="/Applications/Stata/StataSE.app/Contents/MacOS/stata-se")
options("RStata.StataVersion" = 16)

options("RStata.StataPath"="/Applications/Stata/StataBE.app/Contents/MacOS/StataBE -e")
options("RStata.StataVersion" = 18)

command <- "sysuse auto, clear
  sum mpg, d"
stata(command)

Stata('di "Hello World"')
stata('di "Hello World"')
Stata('di 2+2')
stata('clear all')

cmd <- 'import delimited using "Simulations/Products/example.csv"\nreg y x1, nocon\nscpc'

Stata(cmd)
stata("Simulations/Stata/test.do")
Stata("Simulations/Stata/test.do")
## Model testing ------

s_1 <- runif(100) #location vector
cutoff <- quantile(s_1, probs = 0.85)[1]
x_l <- ifelse(s_1<cutoff,-0.15,0.85)
d <- abs(outer(s_1,s_1,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
Sigma_c <- exp(-75.167*d) # Covariance function -> COV matrix
e_l <- mvtnorm::rmvnorm(1,sigma=Sigma_c)
plot(s_1,e_l)
mean(e_l)
mean_pairwise_corr(Sigma_c)

## Searching for c----

obj_f <- function(x,r,d) rho_c(x,d)-r
s_l<-runif(250)
d <- abs(outer(s_l,s_l,"-"))
ans<-uniroot(obj_f,c(0.1,1000),tol=1e-6,c=0.03,d=d)
ans

obj_f <- function(x,r,d) rho_c(x,d)-r

get_c <- function(d){
  ans<-uniroot(obj_f,c(0.1,1000),tol=1e-6,r=0.03,d=d)
  return(ans$root[1])
}

##

load("Simulations/Products/ratio_f_n_chunk_500_.RData")

##


