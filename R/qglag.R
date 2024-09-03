
percentile=function(data,Xnm,q=4){
  for(i in 1:length(Xnm)){
    k=(range(data[,Xnm[i]])[2]-range(data[,Xnm[i]])[1])/q
    kk=seq(range(data[,Xnm[i]])[1],range(data[,Xnm[i]])[2],k)
    data[,Xnm[i]] =(cut(data[,Xnm[i]], breaks = kk, labels = FALSE, include.lowest = TRUE) - 1)
  }
  data
}
quantize<-function (data, expnms, q = 4, breaks = NULL)
{
  e <- new.env()
  e$retbr <- list()
  qt <- function(i) {
    datmat <- as.numeric(unlist(data[, expnms[i]]))
    if (!is.null(breaks)) {
      br <- breaks[[i]]
      e$retbr[[i]] <- breaks[[i]]
    }
    else {
      br <- unique(quantile(datmat, probs = seq(0, 1, by = 1/q),
                            na.rm = TRUE))
      br[1] <- -1e+64
      br[length(br)] <- 1e+64
      e$retbr[[i]] <- br
    }
    cut(datmat, breaks = br, labels = FALSE, include.lowest = TRUE) -
      1
  }
  if (length(expnms) == 1) {
    data[, expnms] <- qt(1)
  }
  else {
    data[, expnms] <- vapply(seq_len(length(expnms)), qt,
                             rep(0, nrow(data)))
  }
  return(list(data = data, breaks = e$retbr))
}


msmfit=function(formula,Xnm,data,cb=NULL,cb_lag=NULL,cb_argvar=NULL,cb_arglag=NULL,
                degree=1){
  #if(is.null(maxlag)){
  #  stop("'maxlag' not recognized")
  #}
  if(is.null(cb)+is.null(cb_lag)+is.null(cb_argvar)+is.null(cb_arglag) !=0&
     is.null(cb)+is.null(cb_lag)+is.null(cb_argvar)+is.null(cb_arglag) !=4){
    stop("cb, cb_lag, cb_argvar, and cb_arglag need to assign values")
  }

  za=as.character(formula)
  zb=str_replace_all(string=strsplit(za[3],split="+",fixed=T)[[1]], pattern=" ", repl="")




  if(!is.null(cb)){
    maxlag=max(cb_lag)
    n1 = length(cb)
    dta=list()
    varr=list()
    lagg=list()
    for(i in 1:n1){
      dta[[i]]=crossbasis(data[,cb[i]], lag = cb_lag[i], argvar = cb_argvar[[i]],
                          arglag = cb_arglag[[i]])
      zb[zb==cb[i]]=paste("dta[[",i,"]]",sep="")
      varr[[i]]=ifelse(identical(attr(dta[[i]],"argvar")$knots,numeric(0)),
                       list(fun=attr(dta[[i]],"argvar")$fun,df=attr(dta[[i]],"df")[1],
                            degree=attr(dta[[i]],"argvar")$degree,Boundary.knots=attr(dta[[i]],"argvar")$Boundary.knots),
                       list(fun=attr(dta[[i]],"argvar")$fun,df=attr(dta[[i]],"df")[1],knots=attr(dta[[i]],"argvar")$knots,
                            degree=attr(dta[[i]],"argvar")$degree,Boundary.knots=attr(dta[[i]],"argvar")$Boundary.knots))

      lagg[[i]]=cb_arglag[[i]]
    }
  }

  dtb = dta


  zc = paste(zb,collapse = "+")
  formula1 = as.formula(paste(za[2],za[1],zc,sep=""))


  model <- glm(formula1,
               family=quasipoisson(),data=data,na.action="na.exclude")

  dt2=NULL
  for(j in 1:100){
    da=data
    da[,Xnm]=(j-1)
    da$psi=(j-1)
    db=da[1:maxlag,]
    db$psi=1000
    da=rbind(db,da)
    #da$y=predict(model,da,type="response")
    dt2=rbind(dt2,da)
  }

  if(!is.null(cb)){
    n1 = length(cb)
    dta=list()
    for(i in 1:n1){
      dta[[i]]=crossbasis(dt2[,cb[i]], lag = cb_lag[i], argvar = varr[[i]],
                          arglag = lagg[[i]])
      colnames(dta[[i]])=paste("dta[[",i,"]]",colnames(dta[[i]]),sep="")
      dt2=cbind(dt2,dta[[i]])
    }
  }



  dt3=cbind(dt2,y=predict(model,dt2,type="response"))
  dt3=dt3[dt3$psi!=1000,]

  uu=glm(y~poly(psi,degree=degree,raw=T),data=dt3,family = quasipoisson())
  names(uu$coefficients)[2:length(names(uu$coefficients))]=paste("psi",1:degree,sep="")
  result=list("ffit"=model,"sfit"=uu,"formula"=formula1,"cb"=dtb,"arg_var"=varr,"lag_var"=lagg)
  result
}




qgdlnm = function(formula,Xnm,data,cb=NULL,cb_lag=NULL,cb_argvar=NULL,cb_arglag=NULL,
                 degree=1,boot=200,sameinterval = TRUE,q=4){


  if(is.null(cb)+is.null(cb_lag)+is.null(cb_argvar)+is.null(cb_arglag) !=0&
     is.null(cb)+is.null(cb_lag)+is.null(cb_argvar)+is.null(cb_arglag) !=4){
    stop("cb, cb_lag, cb_argvar, and cb_arglag need to assign values")
  }

  if(sameinterval==TRUE){
    dt1<-quantize(data, Xnm, q, breaks=NULL)$data
  }else{
    dt1 <- percentile(data, Xnm,q)
  }

  fit=msmfit(formula,Xnm,dt1,cb=cb,cb_lag=cb_lag,cb_argvar=cb_argvar,cb_arglag=cb_arglag,degree=degree)

  sb = as.character(fit$formula)
  formulab = as.formula(paste(sb[2],sb[1],sb[3],sep=""))

  kk=matrix(NA,boot,length(coef(fit$sfit)))
  for(ii in 1:boot){
    sa=sample(1:nrow(dt1),nrow(dt1),replace=T)
    dt11 = dt1[sa,]
    if(!is.null(cb)){
      maxlag=max(cb_lag)
      n1 = length(cb)
      dta=list()
      for(i in 1:n1){
        dta[[i]]=fit$cb[[i]][sa,]
      }
    }

    modela <- glm(formulab,
                 family=quasipoisson(),data=dt11,na.action="na.exclude")

    dt21=NULL
    for(j in 1:100){
      da=dt11
      da[,Xnm]=(j-1)
      da$psi=(j-1)
      db=da[1:maxlag,]
      db$psi=1000
      da=rbind(db,da)
      #da$y=predict(model,da,type="response")
      dt21=rbind(dt21,da)
    }

    if(!is.null(cb)){
      n1 = length(cb)
      dta=list()
      for(i in 1:n1){
        dta[[i]]=crossbasis(dt21[,cb[i]], lag = cb_lag[i], argvar =fit$arg_var[[i]] ,
                            arglag = fit$lag_var[[i]])
        colnames(dta[[i]])=paste("dta[[",i,"]]",colnames(dta[[i]]),sep="")
        dt21=cbind(dt21,dta[[i]])
      }
    }
    dt31=cbind(dt21,y=predict(modela,dt21,type="response"))
    dt31=dt31[dt31$psi!=1000,]

    uu1=glm(y~poly(psi,degree=degree,raw=T),data=dt31,family = quasipoisson())
    names(uu1$coefficients)[2:length(names(uu1$coefficients))]=paste("psi",1:degree,sep="")

    kk[ii,]=uu1$coefficients

  }

  kk=kk[!is.na(kk[,1]),]
  aaa=summary(fit$sfit)
  aaa$cov.scaled[,]=cov(kk)
  aaa$coefficients[,2] = sqrt(diag(cov(kk)))
  aaa$coefficients[,3]=aaa$coefficients[,1]/aaa$coefficients[,2]
  aaa$coefficients[,4]=pnorm(abs(aaa$coefficients[,3]),lower.tail = F)*2

  result = list("summary"=aaa,"fit"=fit$sfit,"coef"=fit$sfit$coefficients,"cov"=cov(kk),"q"=q)
  result
}


###Example
#uzz=qglag(dall ~ anutem +
#        anurh + PM10 + ns(time, 7 * length(unique(year))) +
#        dow,Xnm=c("anutem", "anurh","PM10"),data=dt,cb=c("anutem","anurh"),cb_lag=c(21,5),
#      cb_argvar=list(list(fun=varfun,df=3),list(fun=varfun,df=3)),cb_arglag = list(arglag,arglag),
#      degree=1,boot=20)

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

pg=function(model,cen=50){
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

#pg(uzz,cen=2)

