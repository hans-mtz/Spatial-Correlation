# Loading packages

library(spatstat) #rpoispp
library(fields) #Matern
library(fixest) #feols
# setFixest_nthreads(nthreads=1)  #avoid trouble with parallel
library(mgcv)
# library(broom)
# library(santoku)
library(tidyverse)  
# library(foreach)
# library(parallel)

# library(here)

# source(here("SE_Sims_Code_new.R"))
source("Theta/R-Morgan/SE_Sims_Code_new.R")

# num_cores=detectCores()-2  #number of cores to use

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

basis=8
Structure=0.05
#####Spatial correlation matrix
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


# Morgan's Splines
j=5
sim_y=Sim_y[,j]
sim_x=Sim_x[,j]

gam(sim_y~s(X,bs="bs",k=8, m=2)) |> model.matrix() -> spline_x
gam(sim_y~s(Y,bs="bs",k=8, m=2)) |> model.matrix() -> spline_y
par(mfcol=c(2,1))
plot(rep(X,ncol(spline_x)-1),spline_x[,2:ncol(spline_x)])
plot(rep(Y,ncol(spline_y)-1),spline_y[,2:ncol(spline_y)])

# gam(sim_y~s(X,bs="bs",k=8, m=0)) |> model.matrix() -> spline_x
# gam(sim_y~s(Y,bs="bs",k=8, m=0)) |> model.matrix() -> spline_y
# par(mfcol=c(2,1))
# plot(rep(X,ncol(spline_x)-1),spline_x[,2:ncol(spline_x)])
# plot(rep(Y,ncol(spline_y)-1),spline_y[,2:ncol(spline_y)])


# Morgans HR (fixest)

regs=vector(mode="list", length = 7)
regs[[1]]<-feols(eq_orig,data=df1,vcov="hetero")
regs[[2]]<-feols(eq_orig,data=df1,vcov=vcov_conley(cutoff=cut_off,lat=~Y,lon=~X, vcov_fix = F))
regs[[3]]<-feols(eq_orig,data=df1,vcov=vcov_conley(cutoff=cut_off*2,lat=~Y,lon=~X, vcov_fix = F))
regs[[4]]<-feols(eq_orig,data=df1,vcov=vcov_conley(cutoff=cut_off*3,lat=~Y,lon=~X, vcov_fix = F))
# feols(e=q_orig,data=df1,vcov=vcov_conley(cutoff=500,lat=~I(Y*180),lon=~I(X*360), vcov_fix = F))
regs[[5]]<-feols(eq_orig,data=df1,vcov=vcov_conley(cutoff=5.554,lat=~Y,lon=~X, vcov_fix = F))
regs[[6]]<-feols(eq_orig,data=df1,vcov=vcov_conley(cutoff=5.554*2,lat=~Y,lon=~X, vcov_fix = F))
regs[[7]]<-feols(eq_orig,data=df1,vcov=vcov_conley(cutoff=5.554*3,lat=~Y,lon=~X, vcov_fix = F))

## Closest cutoff 5.5 pixel=0.7
etable(regs)
etable(regs[5:7])


## Get Morgan's Splines

#### Get PCs from basis^2 quadratic b spline


gm_2=bam(sim_y~
            te(X,Y,
              bs="bs",
              k=basis,
              m=2),
          discrete=T)
gm_3=gam(sim_y~
            te(X,Y,
              bs="bs",
              k=basis,
              m=2),
          discrete=T)
df_p=prcomp(model.matrix(gm_2)[,-1])
pc=as.data.frame(df_p$x)     #[,1:best]
spline_tensor_bam<-model.matrix(gm_2)[,-1]
# colnames(spline_tensor_bam)<-paste0("bam.",colnames(spline_tensor_bam))
spline_tensor_gam<-model.matrix(gm_3)[,-1]
# colnames(spline_tensor_gam)<-paste0("gam.",colnames(spline_tensor_gam))


## Morgan's regressions

best=(basis^2)-1  #use all pcs
pc_names=paste("PC",1:best,sep="")   #principal component number to use in each simbest=(basis^2)-1  #use all pcs
eq_orig=as.formula("sim_y~sim_x") 
eq_basis=as.formula(paste("sim_y~sim_x",paste(pc_names,collapse="+"),sep="+")) 
bam_names=paste0("bam",1:best)
colnames(spline_tensor_bam)<-bam_names
gam_names=paste0("gam",1:best)
colnames(spline_tensor_gam)<-gam_names


df1=cbind.data.frame(sim_y,sim_x,pc,X,Y)
df_tst=cbind.data.frame(sim_y,sim_x,pc,X,Y,spline_tensor_bam,spline_tensor_gam)

eq_bam=as.formula(
  paste(
    "sim_y~sim_x+",
     grep("^bam", names(df_tst), value=TRUE) |> paste(collapse="+")
  )
)
eq_bam<-xpd(sim_y~sim_x+regex("bam"),data=df_tst)
eq_gam<-xpd(sim_y~sim_x+regex("gam"),data=df_tst)

## Regressions
morgans_ols_splines<-vector(mode = "list", length = 9)

###HC
morgans_ols_splines[[1]]<-feols(eq_basis,data=df1,
            vcov="hetero")

#Conley SEs with and without spatial basis. Gives 1, 2, 3, and 4 times cutoff chosen. 

morgans_ols_splines[[2]]<-feols(eq_basis,data=df1,
                    vcov=vcov_conley(cutoff=cut_off,lat=~Y,lon=~X, vcov_fix = F))
morgans_ols_splines[[3]]<-feols(eq_basis,data=df1,
                    vcov=vcov_conley(cutoff=2*cut_off,lat=~Y,lon=~X, vcov_fix = F))
morgans_ols_splines[[4]]<-feols(eq_basis,data=df1,
                    vcov=vcov_conley(cutoff=3*cut_off,lat=~Y,lon=~X, vcov_fix = F))
morgans_ols_splines[[5]]<-feols(eq_basis,data=df1,
                    vcov=vcov_conley(cutoff=4*cut_off,lat=~Y,lon=~X, vcov_fix = F))
# The behavior does not change if vcov_fix is true or false, if there is negative values in the variance 
# it gets fixed a la Cameron, Gelbach and Miller (2011)
morgans_ols_splines[[6]]<-feols(eq_bam,data=df_tst,
            vcov="hetero")

morgans_ols_splines[[7]]<-feols(eq_gam,data=df_tst,
            vcov="hetero")

morgans_ols_splines[[8]]<-feols(eq_bam,data=df_tst,
              vcov=vcov_conley(cutoff=cut_off,lat=~Y,lon=~X, vcov_fix = F))

morgans_ols_splines[[9]]<-feols(eq_gam,data=df_tst,
                    vcov=vcov_conley(cutoff=cut_off,lat=~Y,lon=~X, vcov_fix = F))
    

## Results

morgans_ols_splines[seq(2,9, by=2)] #no splines
morgans_ols_splines[seq(3,9, by=2)] #splines
morgans_ols_splines[10:13] #Bam vs Gam
morgans_ols_splines[c(3,12,13)] #Bam vs Gam
# Save data

## To Matlab
write.csv(data.frame(y=sim_y,X=sim_x,s_x=Coords[,1],s_y=Coords[,2]),file="Theta/morgans_data.csv", row.names=FALSE)
write.csv(pc,file="Theta/morgans_splines_pc.csv", row.names=FALSE)
write.csv(spline_tensor_bam,file="Theta/morgans_splines_bam.csv", row.names=FALSE)
write.csv(spline_tensor_gam,file="Theta/morgans_splines_gam.csv", row.names=FALSE)

## To R
save(df1,df_tst,cut_off,Y,X,eq_orig,eq_basis,eq_bam,eq_gam, sim_y, file='Theta/R-Morgan/morgans_data.Rdata')
