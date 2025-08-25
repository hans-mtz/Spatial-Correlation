library(spatstat)
library(fields)
library(fixest)
setFixest_nthreads(nthreads=1)  #avoid trouble with parallel
library(mgcv)
library(broom)
library(santoku)
library(tidyverse)  
library(foreach)
library(parallel)

library(here)

source("Conley_SE.R")
num_cores=detectCores()-2  #number of cores to use

######Correlation between points at distance h apart is Structure*exp(-h/Range). 
#########################
Subjects=500   #number of units
Range=sqrt(2)/10  #Range for Matern correlation on unit square. At dist 0.1 equals Structure/2
Structure=0.75
basis=8

nSim=1000 #number of simulations. Set low at first to check that things are working ok.

##################################
#coordinates of Poisson process
set.seed(123)
pp=rpoispp(2*Subjects)[1:Subjects,]   #points on unit square
#plot(pp)
X=pp$x
Y=pp$y
Coords=cbind(pp$x,pp$y)

#########Spatial Weights for Conley
D = fields::rdist(Coords)
#Wts=ifelse(D<0.05,1,0)                 #rectangular kernel
Wts=dnorm(D,sd=0.05)  /dnorm(0,sd=0.05) #Gaussian weights


#######Simulate noise
KL=Structure*Matern(rdist(x1=Coords),
                    range=Range,
                    smoothness=0.5   #exponential falloff
)+
  diag(nrow(Coords))*(1-Structure) 

#Simulate spatial vars by multiplying cholesky decomp with standard normal variables  
KL=t(chol(KL))  

set.seed(1234)
Sim_y=KL%*%matrix(rnorm(nrow(Coords)*nSim),ncol=nSim)
set.seed(4321)
Sim_x=KL%*%matrix(rnorm(nrow(Coords)*nSim),ncol=nSim)

###set up equations: with and without PCs of spline  
best=(basis^2)-1  #use all pcs
pc_names=paste("PC",1:best,sep="")   #principal component number to use in each sim

eq_orig=as.formula("sim_y~sim_x") 
eq_basis=as.formula(paste("sim_y~sim_x",paste(pc_names,collapse="+"),sep="+")) 



#######simulate SEs for both unadjusted and basis regressions
placebo_sim=function(j){
  sim_y=Sim_y[,j]
  sim_x=Sim_x[,j]
  environment(basis)   
  
  #Get PCs from basis^2 quadratic b spline
  gm_2=bam(sim_y~
             te(X,Y,
                bs="bs",
                k=basis,
                m=2),
           discrete=T)
  df_p=prcomp(model.matrix(gm_2)[,-1])
  pc=as.data.frame(df_p$x)        #[,1:best]
  df1=cbind.data.frame(sim_y,sim_x,pc,X,Y)
  
  c_orig=Conley_SE(lm(eq_orig,data=df1),Wts)[2,4]<0.05
  c_basis=Conley_SE(lm(eq_basis,data=df1),Wts)[2,4]<0.05
  
  cn=data.frame(c_orig,c_basis)
  return(cn)
}

# use regular loop    
#sim_results=list()
# for (i in 1:nSim){
#   sim_results[[i]]=placebo_sim(i)
# }

sim_results=list()
cl <- parallel::makeForkCluster(num_cores)
doParallel::registerDoParallel(cl)
s1=Sys.time()
sim_results=foreach(j=1:nSim) %dopar% {placebo_sim(j)}  
parallel::stopCluster(cl)
Sys.time()-s1
sim_results=plyr::ldply(sim_results)
apply(sim_results,2,mean)
