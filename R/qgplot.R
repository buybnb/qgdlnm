

###Graph

vc=function(coef,cov,k,cen){
  dg=length(coef)
  x=matrix(NA,1,dg)
  for(i in 1:dg){
    x[i]=k^i-cen^i
  }
  sd=sqrt(x%*%cov%*%t(x))
  mm=x%*%coef
  result=list("mean"=mm,"sd"=sd)
  result
}

qgplot=function(model,cen=50){
  rs=matrix(NA,model$q,4)
  for(i in 0:(model$q -1)){
    rs[i+1,1]=i
    kk=vc(coef(model$fit)[2:length(coef(model$fit))],model$cov[2:length(coef(model$fit)),2:length(coef(model$fit))],i,cen=cen)
    rs[i+1,2]=exp(kk$mean)
    rs[i+1,3]=exp(kk$mean - 1.96*kk$sd)
    rs[i+1,4]=exp(kk$mean + 1.96*kk$sd)
  }
  colnames(rs)=c("psi","est","lb","ub")
  rs=as.data.frame(rs)
  gp=ggplot(rs,aes(x=psi,y=est))+geom_line(col="red")+
    xlab("Joint Exposure Quantile")+ylab("Risk Ratio")+theme_bw()+
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.2)
  gp
}

