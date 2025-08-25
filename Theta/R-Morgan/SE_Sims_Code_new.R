performance_se=function(Structure,basis){
  ##simulate fraction of SEs significant at 0.05 for different SE adjustments
  ##basis is number of splines in each direction of tensor 
  
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

##Standard error adjustments below here 
###HC
    HC=feols(eq_orig,data=df1,
                vcov="hetero")
##Small clusters
    small_clust=feols(eq_orig,data=df1,
              vcov=vcov_cluster(cluster=clust_small))
    
#Conley SEs with and without spatial basis. Gives 1, 2, 3, and 4 times cutoff chosen. 
    cnly_1=feols(eq_orig,data=df1,
                 vcov=vcov_conley(cutoff=cut_off,lat=~Y,lon=~X, vcov_fix = F))
    cnly_basis_1=feols(eq_basis,data=df1,
                       vcov=vcov_conley(cutoff=cut_off,lat=~Y,lon=~X, vcov_fix = F))
    cnly_2=feols(eq_orig,data=df1,
                 vcov=vcov_conley(cutoff=2*cut_off,lat=~Y,lon=~X, vcov_fix = F))
    cnly_basis_2=feols(eq_basis,data=df1,
                       vcov=vcov_conley(cutoff=2*cut_off,lat=~Y,lon=~X, vcov_fix = F))
    cnly_3=feols(eq_orig,data=df1,
                 vcov=vcov_conley(cutoff=3*cut_off,lat=~Y,lon=~X, vcov_fix = F))
    cnly_basis_3=feols(eq_basis,data=df1,
                       vcov=vcov_conley(cutoff=3*cut_off,lat=~Y,lon=~X, vcov_fix = F))
    cnly_4=feols(eq_orig,data=df1,
                 vcov=vcov_conley(cutoff=4*cut_off,lat=~Y,lon=~X, vcov_fix = F))
    cnly_basis_4=feols(eq_basis,data=df1,
                       vcov=vcov_conley(cutoff=4*cut_off,lat=~Y,lon=~X, vcov_fix = F))

###BCH with and without basis for specified clusters
#     bch_1=feols(eq_orig,data=df1,
#               vcov=vcov_cluster(cluster=clust_1))
#     bch_basis_1=feols(eq_basis,data=df1,
#               vcov=vcov_cluster(cluster=clust_1))
#     bch_2=feols(eq_orig,data=df1,
#                 vcov=vcov_cluster(cluster=clust_2))
#     bch_basis_2=feols(eq_basis,data=df1,
#                       vcov=vcov_cluster(cluster=clust_2))
#     bch_3=feols(eq_orig,data=df1,
#                 vcov=vcov_cluster(cluster=clust_3))
#     bch_basis_3=feols(eq_basis,data=df1,
#                       vcov=vcov_cluster(cluster=clust_3))
#     bch_4=feols(eq_orig,data=df1,
#                 vcov=vcov_cluster(cluster=clust_4))
#     bch_basis_4=feols(eq_basis,data=df1,
#                       vcov=vcov_cluster(cluster=clust_4))
# 
# ###IM here.
#     im_1=IM_placebo(eq_orig,data=df1,cluster=clust_1)
#     im_basis_1=IM_placebo(eq_basis,data=df1,cluster=clust_1)
#     im_2=IM_placebo(eq_orig,data=df1,cluster=clust_2)
#     im_basis_2=IM_placebo(eq_basis,data=df1,cluster=clust_2)
#     im_3=IM_placebo(eq_orig,data=df1,cluster=clust_3)
#     im_basis_3=IM_placebo(eq_basis,data=df1,cluster=clust_3)
#     im_4=IM_placebo(eq_orig,data=df1,cluster=clust_4)
#     im_basis_4=IM_placebo(eq_basis,data=df1,cluster=clust_4)


 ###Output here. Remember commas if you comment out parts   
    sim_res1=cbind.data.frame(
##Proportion signif at 5%      
#       HC=HC$coeftable[2,4]<0.05 ,
#      clust=small_clust$coeftable[2,4]<0.05 ,
      cnly_1=cnly_1$coeftable[2,4]<0.05 ,
      cnly_basis_1=cnly_basis_1$coeftable[2,4]<0.05 ,
      cnly_2=cnly_2$coeftable[2,4]<0.05 ,
      cnly_basis_2=cnly_basis_2$coeftable[2,4]<0.05 ,
      cnly_3=cnly_3$coeftable[2,4]<0.05 ,
      cnly_basis_3=cnly_basis_3$coeftable[2,4]<0.05 ,
      cnly_4=cnly_4$coeftable[2,4]<0.05 ,
      cnly_basis_4=cnly_basis_4$coeftable[2,4]<0.05 ,
      # bch_1=bch_1$coeftable[2,4]<0.05,
      # bch_basis_1=bch_basis_1$coeftable[2,4]<0.05,
      # bch_2=bch_2$coeftable[2,4]<0.05,
      # bch_basis_2=bch_basis_2$coeftable[2,4]<0.05,
      # bch_3=bch_3$coeftable[2,4]<0.05,
      # bch_basis_3=bch_basis_3$coeftable[2,4]<0.05 ,
      # bch_4=bch_4$coeftable[2,4]<0.05,
      # bch_basis_4=bch_basis_4$coeftable[2,4]<0.05 ,
      # im_1=im_1$p_value<0.05,
      # im_basis_1=im_basis_1$p_value<0.05,
      # im_2=im_2$p_value<0.05,
      # im_basis_2=im_basis_2$p_value<0.05,
      # im_3=im_3$p_value<0.05,
      # im_basis_3=im_basis_3$p_value<0.05,
      # im_4=im_4$p_value<0.05,
      # im_basis_4=im_basis_4$p_value<0.05,
###Now CIs
      ci_HC=4*HC$coeftable[2,2] ,
      ci_clust=4*small_clust$coeftable[2,2] ,
      ci_cnly_1=4*cnly_1$coeftable[2,2] ,
      ci_cnly_basis_1=4*cnly_basis_1$coeftable[2,2] ,
      ci_cnly_2=4*cnly_2$coeftable[2,2] ,
      ci_cnly_basis_2=4*cnly_basis_2$coeftable[2,2] ,
      ci_cnly_3=4*cnly_3$coeftable[2,2] ,
      ci_cnly_basis_3=4*cnly_basis_3$coeftable[2,2] ,
      ci_cnly_4=4*cnly_4$coeftable[2,2] ,
      ci_cnly_basis_4=4*cnly_basis_4$coeftable[2,2]   #,  #add this comma if you uncomment bch and im  
      # ci_bch_1=2*qt(0.975,n_clus[1]-1)*bch_1$coeftable[2,2],
      # ci_bch_basis_1=2*qt(0.975,n_clus[1]-1)*bch_basis_1$coeftable[2,2],
      # ci_bch_2=2*qt(0.975,n_clus[2]-1)*bch_2$coeftable[2,2],
      # ci_bch_basis_2=2*qt(0.975,n_clus[2]-1)*bch_basis_2$coeftable[2,2],
      # ci_bch_3=2*qt(0.975,n_clus[3]-1)*bch_3$coeftable[2,2],
      # ci_bch_basis_3=2*qt(0.975,n_clus[3]-1)*bch_basis_3$coeftable[2,2]  ,
      # ci_bch_4=2*qt(0.975,n_clus[4]-1)*bch_4$coeftable[2,2],
      # ci_bch_basis_4=2*qt(0.975,n_clus[4]-1)*bch_basis_4$coeftable[2,2]  ,
      # ci_im_1=4*im_1$SE,
      # ci_im_basis_1=4*im_basis_1$SE,
      # ci_im_2=4*im_2$SE,
      # ci_im_basis_2=4*im_basis_2$SE,
      # ci_im_3=4*im_3$SE,
      # ci_im_basis_3=4*im_basis_3$SE,
      # ci_im_4=4*im_3$SE,
      # ci_im_basis_4=4*im_basis_3$SE
    )
    return(sim_res1)
  }
  
  cl <- parallel::makeForkCluster(num_cores)
  doParallel::registerDoParallel(cl)
  #s1=Sys.time()
  sim_results=foreach(j=1:nSim) %dopar% {placebo_sim(j)}  
  parallel::stopCluster(cl)
  #Sys.time()-s1
  sim_results=plyr::ldply(sim_results)
  Out=sim_results

  return(Out)
}


IM_placebo=function(eqn,data,cluster){
  #Ibragimov-Mueller estimates
  #t test of null that betas from each cluster are zero
  #pseudo-SE is quarter of confint, pseudo-t is norm variable with same p value
  clust_coefs=split(data,cluster) %>%
    map(~lm(eqn, data = .x)) %>%
    map_df(tidy) %>%
    filter(term == 'sim_x') %>%
    dplyr::select(estimate)
  slopes=as.data.frame(t(as.vector(clust_coefs)))
  names(slopes)=paste("slope_",1:ncol(slopes),sep="")
  if(ncol(slopes)<4){
    slopes$slope_4=NA
  }
   im=t.test(clust_coefs,na.action=na.fail())
  coef_table=data.frame(
    Estimate=im$estimate,
    SE=(im$conf.int[2]-im$conf.int[1])/4,
    t_stat=qnorm(1-im$p.value/2),
    p_value=im$p.value
  )
 
  return(coef_table)
}

