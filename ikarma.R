# Implemented by Fabio M Bayer and Camila Malu da Rosa
# e-mail: bayer@ufsm.br

library(rootSolve)

ikarma<- function (y, ar=NA, ma=NA, link = "logit",diag=1,h=6,X=NA,X_hat=NA,resid=4)
{  
  source("ikarma-mu-phi.R")
  source("ikarma_fit.R")
  
  if (min(y) < 0 || max(y) > 1)
    stop("OUT OF RANGE [0,1]!")
  
  if(is.ts(y)==T)
  {
    freq<-frequency(y)
  }else stop("data can be a time-series object")
  
  
  if(any(is.na(ar))==F) names_phi<-c(paste("phi",ar,sep=""))
  
  if(any(is.na(ma))==F) names_theta<-c(paste("theta",ma,sep=""))
  
  if(any(is.na(X))==F)
  {
    names_beta<-c(paste("beta",1:ncol( as.matrix(X) ),sep=""))
  }
  
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
    stats <- make.link(linktemp)
  else stop(paste(linktemp, "link not available, available links are \"logit\", ",
                  "\"probit\" and \"cloglog\""))
  
  link1 <- structure(list(link = linktemp, 
                          linkfun = stats$linkfun,
                          linkinv = stats$linkinv, 
                          mu.eta = stats$mu.eta, 
                          diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  )
  )
  
  fit1 <- ikarma.fit(y, ar, ma, link1, names_phi, names_theta, names_beta, diag, h, X, X_hat,resid=resid) # model estimation
  
  return(fit1)
}


