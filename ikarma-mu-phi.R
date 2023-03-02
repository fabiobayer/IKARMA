# Implemented by Fabio M Bayer and Camila Malu da Rosa
# e-mail: bayer@ufsm.br

library(extraDistr)

dikarma<-function(y,lambda,b,mu,prec) 
{
  f0 = ifelse((y == 0), lambda*(1-b), 0)
  f1 = ifelse((y == 1), lambda*b, 0)
  
  f =  ifelse((y == 0), 0, (1-lambda)*dkumar(y,prec,log(0.5)/log(1-mu^prec))) 
  f =  ifelse((y == 1), 0, f)
  
  (f0 + f1 + f)
}

pikarma<-function(y,lambda,b,mu,prec)
{
  b0 = lambda*(1-b)
  b1 = ifelse((y == 1), lambda*(b), 0)
  b<- (1-lambda)*(1-(1-(y^prec))^( (log(0.5))/(log(1-(mu^prec))) ))
  ret = ifelse(y==1, b0+b1+b, b0+b)
  ret
}

qikarma<-function(u,lambda, b, mu,prec)
{
  n = length(u)
  if(length(lambda)==1) lambda <- rep(lambda,n)
  if(length(b)==1) b <- rep(b,n)
  if(length(mu)==1) mu <- rep(mu,n)
  if(length(prec)==1) prec <- rep(prec,n)
  if(n>1)
  {
    ret=rep(NA,n)
    for(i in 1:n)
    {
      if(u[i] < (lambda[i]*(1-b[i])) ) 
      {
        ret[i]=0
      }else{
        if(u[i] > (1-lambda[i]*b[i]) ) 
        {
          ret[i]=1
        }else{
          ret[i]= qkumar((u[i]-lambda[i]*(1-b[i]))/(1-lambda[i]), prec[i], log(0.5)/log(1-mu[i]^prec[i]))
        }         
      }
    }
  }else{
    if(u < (lambda*(1-b)) ) 
    {
      ret=0
    }else{
      if(u > (1-lambda*b) ) 
      {
        ret=1
      }else{
        ret=qkumar((u-lambda*(1-b))/(1-lambda), prec, log(0.5)/log(1-mu^prec)) 
      }     
    }
  }
  ret
}

rikarma<-function(n,lambda,b,mu,prec)
{
  u <- runif(n)
  qikarma(u,lambda=lambda,b=b,mu=mu,prec=prec)
}

meus <- function(yt,b=NA,n=NA,link)
{
  nyt <- length(yt)
  linkfun <- link$linkfun
  
  s <- abs(yt[nyt]-(1-b))
  sstar <- pmin(pmax(s, 0.5/n), (n-0.5)/n)
  snorm <- linkfun(sstar)
  
  return(snorm)
}
