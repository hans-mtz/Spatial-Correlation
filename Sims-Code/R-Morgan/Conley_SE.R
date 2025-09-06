
########conley standard errors
Conley_SE=function(ols,W){
###ols is lm object, W is weighting matrix
  X=qr.X(ols$qr)
  N=nrow(X)
  k=ncol(X)
  U=diag(ols$res)         
  V=t(X)%*%(U%*%W%*%U)%*%X/N
  xx.inv=N*chol2inv(qr.R(ols$qr))
  Cov=xx.inv%*%V%*%xx.inv/N
  conley.se=sqrt(diag(Cov))
  conley.t=summary(ols)$coef[,1]/conley.se
  conley.p=2*(1-pnorm(abs(conley.t)))
  conley=cbind.data.frame(coef=ols$coef,conley.se,conley.t,conley.p)
  return(conley)  #adjusted standard errors
}
