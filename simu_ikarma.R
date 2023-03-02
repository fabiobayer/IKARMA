# Implemented by Fabio M Bayer and Camila Malu da Rosa
# e-mail: bayer@ufsm.br

# n: sample size
# freq: anual frequency

library(gamlss)
library(extraDistr)

simu.ikarma <- function(n, alpha=0.5,beta=NA,
                        phi=NA,theta=NA, prec=10,
                        omega1, omega2, b,
                        X=NA, freq=12,link="logit")
{
  source("ikarma-mu-phi.R")
  
  ar<-NA
  ma<-NA
  
  if(any(is.na(phi)==F))
  {
    ar <- 1:length(phi)
  }
  
  if(any(is.na(theta)==F))
  {
    ma <- 1:length(theta)
  }
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
  {  
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"logit\", ",
               "\"probit\" and \"cloglog\""))
  } 
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv
  )
  )
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
  ###### IKARMA model sem regressores para (0,1] 

  if(any(is.na(phi)==F) && any(is.na(theta)==F) && any(is.na(X)==T) && b==1) #n0==0 && n1!=0)
  {
    
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)
    
    ynew <-rep(alpha,(n+m))
    gama <- linkinv(ynew)
    
    
    error<-rep(0,n+m) # E(error)=0 
    eta1<- eta2 <- s <- lambda <-mu <- rep(0,n+m)
    y <- y1 <- ys <- rep(0.5,n+m)
    
    for(i in (m+1):(n+m))
    {
      s[i-1] <- meus(y[1:(i-1)],b,n,link)
      
      eta1[i]<- omega1 + omega2*s[i-1]
      eta2[i]  <- alpha + as.numeric(phi%*%linkfun(ys[i-ar])) + as.numeric(theta%*%error[i-ma])
      
      lambda[i] <- linkinv(eta1[i])
      mu[i] <- linkinv(eta2[i])
      y[i]    <- rikarma(1,lambda[i],b,mu[i],prec) 
      ys[i]<- pmin(pmax(y[i], 0.5/n), (n-0.5)/n)
      error[i] <- linkfun(ys[i])-eta2[i]
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } 
  
  ##################################################################################3
  ###### IKARMA model sem regressores para [0,1) 
  
  if(any(is.na(phi)==F) && any(is.na(theta)==F) && any(is.na(X)==T) && b==0) 
    
  {
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)
    
    ynew <-rep(alpha,(n+m))
    gama <- linkinv(ynew)
    
    
    error<-rep(0,n+m) 
    eta1<- eta2 <- s <- lambda <-mu <- rep(0,n+m)
    y <-y1 <- ys <- rep(0.5,n+m)
   
    for(i in (m+1):(n+m))
    {
      s[i-1] <- meus(y[1:(i-1)],b,n,link)
      eta1[i]<- omega1 + omega2*s[i-1]
      eta2[i]  <- alpha + as.numeric(phi%*%linkfun(ys[i-ar])) + as.numeric(theta%*%error[i-ma])
      
      lambda[i] <- linkinv(eta1[i])
      mu[i] <- linkinv(eta2[i])
      y[i]    <- rikarma(1,lambda[i],b,mu[i],prec)
      ys[i]<- pmin(pmax(y[i], 0.5/n), (n-0.5)/n)
      error[i] <- linkfun(ys[i])-eta2[i]
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  }
  
  
  # ############################################################################################
  # ###### IKARMAX model com regressores e (0,1]
  
  if(any(is.na(phi)==F) && any(is.na(theta)==F) && any(is.na(X)==F) && b==1)
  {
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)

    ynew <-rep(alpha,(n+m))
    gama <- linkinv(ynew)

    error<-rep(0,n+m) 
    eta1<- eta2 <- s <- lambda <-mu <- rep(0,n+m)
    y <-y1 <- ys <- rep(0.5,n+m)
    
    X<-as.matrix(X)
    for(i in (m+1):(n+m))
    {
      s[i-1] <- meus(y[1:(i-1)],b,n,link)
      eta1[i]<- omega1 + omega2*s[i-1]
      eta2[i]  <- alpha + as.numeric(X[(i),]%*%as.matrix(beta)) + 
        as.numeric(phi%*%(linkfun(ys[i-ar])- X[(i-ar),]%*%as.matrix(beta))) + as.numeric(theta%*%error[i-ma])
      lambda[i] <- linkinv(eta1[i])
      mu[i] <- linkinv(eta2[i])
      y[i]    <- rikarma(1,lambda[i],b,mu[i],prec) 
      ys[i]<- pmin(pmax(y[i], 0.5/n), (n-0.5)/n)
      error[i] <- linkfun(ys[i])-eta2[i]
    }
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  }
  
  # ############################################################################################
  # ###### IKARMAX model com regressores e [0,1)
  
  if(any(is.na(phi)==F) && any(is.na(theta)==F) && any(is.na(X)==F) && b==0)
  {
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)
    
    ynew <-rep(alpha,(n+m))
    gama <- linkinv(ynew)
    
    error<-rep(0,n+m) 
    eta1<- eta2 <- s <- lambda <-mu <- rep(0,n+m)
    y <-y1 <- ys <- rep(0.5,n+m)
    
    X<-as.matrix(X)
    for(i in (m+1):(n+m))
    {
      s[i-1] <- meus(y[1:(i-1)],b,n,link)
      eta1[i]<- omega1 + omega2*s[i-1]
      eta2[i]  <- alpha + as.numeric(X[(i),]%*%as.matrix(beta)) + 
        as.numeric(phi%*%(linkfun(ys[i-ar])- X[(i-ar),]%*%as.matrix(beta))) + as.numeric(theta%*%error[i-ma]) ###FB Havia parenteses no local errado
      lambda[i] <- linkinv(eta1[i])
      mu[i] <- linkinv(eta2[i])
      y[i]    <- rikarma(1,lambda[i],b,mu[i],prec) 
      ys[i]<- pmin(pmax(y[i], 0.5/n), (n-0.5)/n)
      error[i] <- linkfun(ys[i])-eta2[i]
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  }
 
}

