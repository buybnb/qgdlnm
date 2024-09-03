\name{qgdlnm}
\alias{qgdlnm}
\title{Quantile G-computation Distributed Lag Non-linear Model (QG-DLNM)}
\description{
A method that incorporates the distributed lag non-linear model (DLNM) into quantile g-computation, allowing it to effectively manage lag effects.
}
\usage{
qgdlnm(formula,Xnm,data,cb=NULL,cb_lag=NULL,cb_argvar=NULL,cb_arglag=NULL,
                 degree=1,boot=200,sameinterval = TRUE,q=4)
}
\arguments{
  \item{formula}{R style formula: a symbolic description of the model to be fitted}
  \item{data}{ an input dataframe}
  \item{Xnm}{ character vector of exposures of interest}
  \item{cb}{character vector of exposures of interest with crossbasis}
  \item{cb_lag}{maximum lag for crossbasis exposure}
  \item{cb_argvar,cb_arglag}{ lists of arguments to be passed to the function onebasis for generating the two basis matrices for predictor and lags, respectively.}
  \item{degree}{ polynomial bases for marginal model (e.g. degree = 2 allows that the relationship between the whole exposure mixture and the outcome is quadratic (default = 1).}
  \item{boot}{integer: number of bootstrap iterations}
  \item{q}{Number of quantiles used to create quantile indicator variables representing the exposure variables.}
  \item{sameinterval}{quantiles with same range or not (default = T).}

}

\value{
  An object of class \code{qglag}, which is a list with the following components:
  \item{summary,fit}{Summary of model output}
  \item{coef}{A vector of coefficients}
  \item{cov}{A vector of covariance.}
  \item{q}{Number of quantiles}
}
\details{
  This code include Breimans' unsupervised method and Siegel and her colleagues' purposeful clustering method to calculate initial labels
  To imput user specified initial labels, please use the function \code{initial}
}
\references{
  Lin Z, Wang M, Ma J, et al. The joint effects of mixture exposure to multiple meteorological factors on step count: A panel study in China[J]. \emph{Environmental Pollution}, 2024, 346: 123469.
}
\examples{

dd=chicagoNMMAPS[1:100,]
varfun = "bs"
arglag <- list(fun="bs",df=4)
uzz=qgdlnm(death ~ temp +
            rhum + pm10 + o3 + ns(time, 7 * length(unique(year))) +
        dow,Xnm=c("temp", "rhum","pm10"),data=dd,cb=c("temp","rhum"),cb_lag=c(21,5),
      cb_argvar=list(list(fun=varfun,df=3),list(fun=varfun,df=3)),cb_arglag = list(arglag,arglag),
      degree=1,boot=20)
uzz

pg(uzz,cen=2)

}
