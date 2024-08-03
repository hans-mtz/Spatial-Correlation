library(spatstat)
library(fields)
library(tidyverse)  
library(foreach)
library(parallel)
library(abind)

library(here)

source(here("Conley_SE.R"))

num_cores=detectCores()-2  #number of processors to use


#########################
Subjects=500   #number of units
Range=sqrt(2)/10
Structure=0.75

spline=8 

nSim=100

##################################
#coordinates
set.seed(123)
pp=rpoispp(2*Subjects)[1:Subjects,]   #points on unit square

X=pp$x
Y=pp$y
Coords=cbind(pp$x,pp$y)

#########Spatial Weights for Conley
D = fields::rdist(Coords)               #distance between points on unit square
#Wts=ifelse(D<0.05,1,0)                 #rectangular kernel
Wts=dnorm(D,sd=0.05)  /dnorm(0,sd=0.05) #Gaussian weights

#####nearest neighbours for residuals
nearest=knn2nb(knearneigh(Coords,k=1,longlat = F))   #k nearest
nn=unlist(nearest)      #nearest neighbour of each point
    
########Simulate outcome
#clust_4=factor(cluster::pam(Coords,4)$cluster)
#performance_se(0.5)

# CW = rdist(Coords)  #Conley weights Gaussian
# CW=dnorm(CW,sd=0.05)/dnorm(0,sd=0.05)

#DistMat=dist_mat(data.frame(X=X,Y=Y), lat = "Y", lon = "X")   #distances for conleyreg



##simulate fraction of SEs significant at 0.05 for different SE adjustments
KL=Structure*Matern(rdist(x1=Coords),
                    range=Range,
                    smoothness=0.5   #exponential falloff
)+
  diag(nrow(Coords))*(1-Structure) 

KL=t(chol(KL))  

set.seed(1234)
Sim_y=KL%*%matrix(rnorm(nrow(Coords)*nSim),ncol=nSim)
set.seed(4321)
Sim_x=KL%*%matrix(rnorm(nrow(Coords)*nSim),ncol=nSim)

spline_sim=function(j){
  sim_y=Sim_y[,j]
  sim_x=Sim_x[,j]
  
  gm_2=bam(sim_y~
             te(X,Y,bs="bs",
                k=spline,
                m=1),
           discrete=T)
  df_p=prcomp(model.matrix(gm_2))
  df_p=as.data.frame(df_p$x)
  
  conley_test=list()
  for (k in 2:ncol(df_p)){
    pc=df_p[,1:k]
    pc_names=paste("PC",1:ncol(pc),sep="")   #principal component number to use in each sim
    eq_basis=as.formula(paste("sim_y~sim_x",paste(pc_names,collapse="+"),sep="+")) 
    eq_dep=as.formula(paste("sim_y~",paste(pc_names,collapse="+"),sep="+")) 
    df1=cbind.data.frame(sim_y,sim_x,pc,X,Y)
    BIC=BIC(lm(eq_dep,data=df1))
    unadj=summary(lm(eq_basis,data=df1))
    cnly_basis=Conley_SE(lm(eq_basis,data=df1),Wts)
    sim_res1=cbind.data.frame(
      BIC,
      unadj_p=unadj$coefficients[2,4]<0.1 ,
      cnly_basis=cnly_basis[2,4]<0.1,
      cnly_width=4*cnly_basis[2,2],
      nearest_cor=cor(cbind.data.frame(unadj$residuals,unadj$residuals[nn]))[2,1]
    )
    conley_test[[k]]=sim_res1
  }
  conley_test=list_rbind(conley_test)
  
} 

# s1=Sys.time()   #check it is working and find out how long it takes
# aa=spline_sim(3)
# Sys.time()-s1

sim_results=list()
cl <- parallel::makeForkCluster(num_cores)
doParallel::registerDoParallel(cl)
s1=Sys.time()
sim_results=foreach(j=1:nSim) %dopar% {spline_sim(j)}  
parallel::stopCluster(cl)
Sys.time()-s1


all.matrix <- abind(sim_results, along=3)
hac_basis= data.frame(apply(all.matrix, c(1,2), mean, na.rm=T)  )
hac_basis$index=2:(spline^2)

hac_basis[1:nrow(hac_basis)-1,] %>%    ###Remove last row as X'X here are mostly (all?) singular
  select(
    index,
    BIC,
    "Het Reject"=unadj_p,
         "HAC 0.1 Reject"=cnly_basis,
         "Avg. CI"=cnly_width,
         "Neighbor Corr"=nearest_cor
  ) %>% 
  #select(-unadj_p) %>% 
  pivot_longer(-index) %>% 
  mutate(name=factor(name,levels=c("HAC 0.1 Reject","Avg. CI","Het Reject","BIC" ,"Neighbor Corr")))%>% 
  ggplot(aes(index,value)) +
  geom_line()+
  facet_wrap(vars(name),scales="free",nrow=2)+
  theme_bw(base_size = 15)+
  labs(subtitle=paste("Properties of Alternate G Specifications Using Principal Components."),
       y="",x="Basis dimension.")


#ggsave(here("HAC_triangular_kernel.pdf"),width=8,height=6)




