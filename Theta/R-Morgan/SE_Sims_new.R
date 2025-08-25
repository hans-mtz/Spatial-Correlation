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

source(here("SE_Sims_Code_new.R"))
source("Theta/R-Morgan/SE_Sims_Code_new.R")
num_cores=detectCores()-2  #number of cores to use

######Correlation between points at distance h apart is Structure*exp(-h/Range). 
#########################
Subjects=500   #number of units
Range=sqrt(2)/10  #Range for Matern correlation on unit square. At dist 0.1 equals Structure/2

#cut_off=c(5,10,15,20)*1.1   #conley cutoff in kilometres (unit square is in longitude and latitude: 110 km so multiply by 1.1).
cut_off=5.5   #in km. Equals 0.05 of square of latitude equal to 110km
nSim=10   #number of simulations. Set low at first to check that things are working ok.

##################################
#coordinates of Poisson process
set.seed(123)
pp=rpoispp(2*Subjects)[1:Subjects,]   #points on unit square
#plot(pp)
X=pp$x
Y=pp$y
Coords=cbind(pp$x,pp$y)

#######k-medoids clusters for BCH and IM. Change number as desired
n_clus=c(4,6,8,10)     #desired clusters
clust_1=factor(cluster::pam(Coords,n_clus[1])$cluster)
clust_2=factor(cluster::pam(Coords,n_clus[2])$cluster)
clust_3=factor(cluster::pam(Coords,n_clus[3])$cluster)
clust_4=factor(cluster::pam(Coords,n_clus[4])$cluster)

######square grid of small clusters for regular clustered SEs
sm_clus=5  #number each direction
cuts_y=santoku::chop_quantiles(Coords[,1], (1:(sm_clus-1))/sm_clus) 
cuts_x=santoku::chop_quantiles(Coords[1,], (1:(sm_clus-1))/sm_clus) 
clust_small=factor(paste(cuts_x,cuts_y,sep=""))


#########Simulate outcome for different values of structure rho and different basis numbers
#aa=performance_se(Structure=0.5,basis=8)  #perform quick tryout to make sure it is working
s1=Sys.time()
out_sim=list()
rho=seq(0,1,by=0.1)   #Structures to try
for (i in 1:length(rho)){
  out_sim[[i]]=performance_se(Structure=rho[i],basis=8)
}
Sys.time()-s1
#save(out_sim,file=here("results","all_sims.Rda"))
#######################Average over simulations to get rejection rates and average CIs
avges=list()
for (j in 1:length(rho)){
avges[[j]]=data.frame(t(apply(out_sim[[j]],2,mean)))
}

tabs=list_rbind(avges) %>%    
  mutate(corr=rho/2) %>% #add correlation at distance 0.1
  relocate(corr)
View(tabs)
#save(tabs_long_range_adj,file=here("results","long_range_adj_cnly_sims.Rda"))

#######Now make tables for CIs and rejection rates, with and without splines
#######Select cnly bch im etc
#Tables of rejection rates at 5%
no_basis_signif=tabs %>%     
  select(-contains("_basis")) %>% select(-contains("ci")) %>% 
  select(corr,contains("cnly"))
basis_signif=tabs %>% 
  select(corr,contains("_basis"),-contains("ci"))%>% 
  select(corr,contains("cnly"))
signif=left_join(no_basis_signif,basis_signif,by="corr")

#Tables of CIs
no_basis_ci=tabs %>%      #table of CIs
  select(corr,contains("ci"),-contains("_basis")) %>% 
  select(corr,contains("cnly"))
basis_ci=tabs %>%      #table of CIs
  select(corr,contains("ci")) %>% select(contains("corr"),contains("_basis"))%>% 
  select(corr,contains("cnly"))
ci=left_join(no_basis_ci,basis_ci,by="corr")

print(xtable::xtable(signif,digits=2),include.rownames=FALSE,)  #latex output
print(xtable::xtable(ci,digits=2),include.rownames=FALSE,)  #latex output



####################################plot results
