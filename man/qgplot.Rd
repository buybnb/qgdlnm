\name{qgplot}
\alias{qgplot}
\title{Quantile G-computation Distributed Lag Non-linear Model Plot}
\description{
Draw the mixture effect plot from qgdlnm() function
}
\usage{
qgplot(model,cen=50)
}
\arguments{
  \item{model}{model output from qgdlnm() function}
  \item{cen}{The reference value used as the counterfactual scenario (the comparator)}


}


\details{
  qgplot() is used to create mixture effect plot from qgdlnm() function
}
\references{
  Lin Z, Wang M, Ma J, et al. The joint effects of mixture exposure to multiple meteorological factors on step count: A panel study in China[J]. \emph{Environmental Pollution}, 2024, 346: 123469.
}
\examples{
library(dlnm)
dd=chicagoNMMAPS[1:100,]
varfun = "bs"
arglag <- list(fun="bs",df=4)
uzz=qgdlnm(death ~ temp +
            rhum + pm10 + o3 + ns(time, 7 * length(unique(year))) +
        dow,Xnm=c("temp", "rhum","pm10"),data=dd,cb=c("temp","rhum"),cb_lag=c(21,5),
      cb_argvar=list(list(fun=varfun,df=3),list(fun=varfun,df=3)),cb_arglag = list(arglag,arglag),
      degree=1,boot=20)
uzz

qgplot(uzz,cen=2)

}
