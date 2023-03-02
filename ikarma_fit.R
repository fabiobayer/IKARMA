# Implemented by Fabio M Bayer and Camila Malu da Rosa
# e-mail: bayer@ufsm.br

ikarma.fit<- function (y, ar, ma, link, names_phi,names_theta,names_beta,diag,h1,X,X_hat,resid)
{
  maxit1<-1000
  
  n<- length(y) # sample size
  n0<- sum(y==0) # number of zeros
  n1<- sum(y==1) # number of ones
  
  y_bar<- mean(y)
  
  delta0<- n0/n
  delta1<- n1/n
  
  a0<- delta0/(1-y_bar)
  a1<- delta1/(y_bar)
  
  ys01<- y
  if(any(y==0)) ys01<- ys01[-which(y==0)]
  if(any(y==1)) ys01<- ys01[-which(ys01==1)]
  
  y01<-y
  y01[which(y01==0)]<- min(ys01)
  y01[which(y01==1)]<- max(ys01)
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta <-  link$mu.eta
  diflink <- link$diflink
  
  ynew01 = linkfun(y01)
  ystar01 = log(y01/(1-y01))
  
  p <- max(ar)
  q <- max(ma)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)
  
  y1<-y[(m+1):n]
  ys = pmin(pmax(y, 0.5/n), (n-0.5)/n)
 
  if(any(is.na(ar)==F)) 
  {
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    
    for(i in 1:(n-m))
    {
      P[i,] <- linkfun(ys[i+m-ar])
    }
    
    Z <- cbind(rep(1,(n-m)),P)
  }else{
    Z <- as.matrix(rep(1,(n-m)))
  }
 
  if(any(is.na(X)==T)) 
  {
    x <- as.matrix(Z)
    Y01 <- y01[(m+1):n]
    Ynew01 = linkfun(Y01)
    Ystar01 = log(Y01/(1-Y01))
    ajuste = lm.fit(x, Ynew01)
    mqo = c(ajuste$coef)
    k = length(mqo)
    ntemp = length(Y01)
    mean = fitted(ajuste)
    mean = linkinv(mean)
    dlink = diflink(mean)
    er = residuals(ajuste)
    sigma2 = sum(er^2)/((ntemp - k) * (dlink)^2)
    prec = 1/ntemp * sum(mean * (1 - mean)/sigma2 - 1)
  }else{ 
    X_hat<-as.matrix(X_hat)
    X<-as.matrix(X)
    x <- cbind(as.matrix(Z),X[(m+1):n,])
    Y01 <- y01[(m+1):n]
    Ynew01 = linkfun(Y01)
    Ystar01 = log(Y01/(1-Y01))
    ajuste = lm.fit(x, Ynew01)
    mqo = c(ajuste$coef)
    k = length(mqo)
    ntemp = length(Y01)
    mean = fitted(ajuste)
    mean = linkinv(mean)
    dlink = diflink(mean)
    er = residuals(ajuste)
    sigma2 = sum(er^2)/((ntemp - k) * (dlink)^2)
    prec = 1/ntemp * sum(mean * (1 - mean)/sigma2 - 1)
    
  }

  z <- c()
  
  ###############################################################################
  # ######### ARMA model sem regressores para (0,1] - IKARMA
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==T) && n0==0 && n1!=0)
  { 
    b <- 1
    reg <- c(linkfun(y_bar),rep(0,p1+q1), 0.2*prec, linkfun(delta1), 0) 
    
    loglik <- function(d) 
    {
      alpha <- d[1]
      phi = d[2:(p1+1)] 
      theta = d[(p1+2):(p1+q1+1)]
      prec <- d[p1+q1+2] 
      omega1 <- d[p1+q1+3]
      omega2 <- d[p1+q1+4]
      
      error<-rep(0,n) 
      eta1<- eta2<- s <- rep(NA,n)
      
      for(i in (m+1):n)
      {
        s[i-1] <- meus(yt=y[1:(i-1)], b=b, n=n,link)
        eta1[i]<- omega1 + omega2*s[i-1]
        eta2[i] <- alpha + (phi%*%linkfun(ys[i-ar])) + (theta%*%error[i-ma])
        error[i] <- linkfun(ys[i])-eta2[i]
      }
      lambda <- linkinv(eta1[(m+1):n]) 
      mu <- linkinv(eta2[(m+1):n]) 
      
      l1 <- ifelse( (y1==0) | (y1==1), log(lambda),0) +
            ifelse( (y1==0) | (y1==1), 0,log(1-lambda)) 
      
      l2 <- ifelse((y1==0) |(y1==1), 0, 
                   log(prec) + 
                     ( log(log(0.5)/log(1-mu^prec)) ) +
                     (prec-1)*(log(y1)) +
                     ( (log(0.5)/log(1-mu^prec))-1)* ( log(1-y1^prec) ) )
      
      l3a <- ifelse( (y1==1),log(b),0)
      l3b <- ifelse( (y1==0),log(1-b),0)
      l3 <- l3a + l3b
    
      sum(l1 + l2+ l3)
    } 
    
    escore <- function(w)
    {
      alpha <- w[1]
      phi = w[2:(p1+1)]
      theta = w[(p1+2):(p1+q1+1)]
      prec <- w[p1+q1+2] 
      omega1 <- w[p1+q1+3]
      omega2 <- w[p1+q1+4]
      
      error<-rep(0,n) 
      eta1<- eta2<- s <- rep(NA,n)

      for(i in (m+1):n)
      {
        s[i-1] <- meus(yt=y[1:(i-1)], b=b, n=n,link)
        eta1[i]<- omega1 + omega2*s[i-1]
        eta2[i] <- alpha + (phi%*%linkfun(ys[i-ar])) + (theta%*%error[i-ma])
        error[i] <- linkfun(ys[i])-eta2[i]
      }
      
      lambda <- linkinv(eta1[(m+1):n])
      mu <- linkinv(eta2[(m+1):n])
      y1<-y[(m+1):n]
      s <- s[m:(n-1)]

      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        R[i,] <- error[i+m-ma]
      }
      
      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      }
      
      a <- ifelse( (y1 == 0) , (1)/(lambda), (-1)/((1-lambda)))
      a <- ifelse( (y1 == 1) , (1)/(lambda), a)
      
      I <- as.matrix(rep(1,n-m))
      
      delta <- log(0.5)/log(1-mu^prec)
      
      ci <- ifelse((y1 == 0) | (y1 == 1), 0, 
                   ( (mu^(prec-1)) / ((1-mu^prec)*log(1-mu^prec)) )*( delta*log(1-y1^prec)+1 ))
      c <- prec*ci
      
      nu <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rR <- deta.dtheta[(m+1):n,]
      
      mT1 <- diag(mu.eta(eta1[(m+1):n]))
      mT2 <- diag(mu.eta(eta2[(m+1):n]))
      
      phi1 <- 1/prec
      phi2 <- log(y1)
      phi3 <- ci*mu*log(mu)
      phi4 <- delta-1
      phi5 <- ( (y1^prec)*log(y1) ) / (1-y1^prec)
      dphi <- as.vector(phi1+phi2+phi3-(phi4*phi5)) 
      dphi1 <- ifelse( (y1==0)|(y1==1), 0,dphi) 
      
      Ualpha <- t(nu) %*% mT2 %*% c
      Uphi <-   t(rP) %*% mT2 %*% c
      Utheta <- t(rR) %*% mT2 %*% c
      Uprec <- sum(dphi1)
      Uomega1 <- t(a) %*% mT1 %*% I
      Uomega2 <- t(a) %*% mT1 %*% s

      rval <- c(Ualpha, Uphi, Utheta, Uprec, Uomega1, Uomega2)
    }
    names_par <- c("alpha", names_phi, names_theta, "precision","omega1", "omega2")
    
    opt <- optim(reg, loglik, escore, 
                 method = "BFGS", 
                 control = list(fnscale = -1, reltol = 1e-12,maxit=200))

    if (opt$conv != 0) 
      warning("FUNCTION DID NOT CONVERGE!")
    
    z$conv <- opt$conv
    coef <- (opt$par)
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    prec <- coef[p1+q1+2] 
    omega1 <- coef[p1+q1+3]
    omega2 <- coef[p1+q1+4]
    
    z$alpha<- alpha
    z$phi <- phi
    z$theta <- theta
    z$prec <- prec
    z$omega1 <- omega1 
    z$omega2 <- omega2 
   
    errorhat<-rep(0,n) 
    etahat1 <- etahat2 <- s <- rep(NA,n)
    
    for(i in (m+1):n)
    {
      s[i-1] <- meus(yt=y[1:(i-1)], b=b, n=n,link)
      etahat1[i] <- omega1 + omega2*s[i-1]
      etahat2[i] <-alpha + (phi%*%linkfun(ys[i-ar])) + (theta%*%errorhat[i-ma])
      errorhat[i] <- linkfun(ys[i])-etahat2[i]
    }
    
    lambdahat <- linkinv(etahat1[(m+1):n])
    muhat <- linkinv(etahat2[(m+1):n])
    y1<-y[(m+1):n]
    s <- s[m:(n-1)]
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat1 <- etahat1
    z$etahat2 <- etahat2
    z$errorhat <- errorhat
    z$s <- s
    
    z$mu <- muhat
    z$lambda <- lambdahat
    z$pzero <- lambdahat*(1-b)
    z$pum <- lambdahat*(b)

    Q <- cbind(1,s)
   
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      R[i,] <- errorhat[i+m-ma]
    }
    
    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    
    for(i in (m+1):n)
    {
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
    }
    
    L <- diag(-1/((lambdahat)*(1-lambdahat)))
   
    v <- matrix(as.vector(deta.dalpha[(m+1):n]),ncol=1)
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    
    mT1 <- diag(mu.eta(etahat1[(m+1):n]))
    mT2 <- diag(mu.eta(etahat2[(m+1):n]))
  
    delta <- log(0.5)/log(1-muhat^prec)
    
    ci <- ifelse((y1 == 0) | (y1 == 1), 0, 
                 ((muhat^(prec-1))/((1-muhat^prec)*log(1-muhat^prec)))*(delta*log(1-y1^prec)+1))
    
    Vc <- prec*ci
    
    zeta1 <- (muhat^(prec-2))/((1-muhat^prec) * (log(1-muhat^prec)))
    zeta2 <- (muhat^(2*prec-2))/((1-muhat^prec)^2 * (log(1-muhat^prec))^2)
    
    W <- diag(as.vector((lambdahat -1)*  prec^2 * zeta2))
    vI <- matrix(rep(1,(n-m)),ncol=1)
    euler <- 0.5772156649
    
    d <- -prec*muhat*log(muhat)*zeta2-prec*muhat*delta*zeta1*((1-digamma(delta+1)-euler)/((delta-1)*prec))
    D <- diag(as.vector(d))
    
    N1 <- 1/(prec^2) 
    N2<- delta*(muhat^2)*zeta2*(log(muhat)^2)*log(1-y1^prec)
    N3<- 2*delta*(muhat^2)*log(muhat)*zeta1*(( (y1^prec) * log(y1))/(1-y1^prec))
    N4<- (delta-1)*( ((y1^prec)*log(y1)^2)/((1-y1^prec)^2) )   
    N5<- ci*muhat*(log(muhat)^2)*( (zeta1*(muhat^2)) +(1/(1-muhat^prec)) )
    
    N<- ifelse((y1 == 0), 0,-N1+N2-N3-N4+N5)
    N<- ifelse((y1 == 1), 0, N)
    
    o = t(matrix(data = 0, nrow = 1, ncol = 2))
    op = t(matrix(data = 0, nrow = p1, ncol = 2))
    ot = t(matrix(data = 0, nrow = q1, ncol = 2))
    
    Koo <- t(Q) %*% L %*% mT1^2 %*% Q
    
    Kaa <- t(v) %*% W %*% mT2^2 %*% v
    Kap <- t(v) %*% W %*% mT2^2 %*% rP
    Kpa <- t(Kap)
    Kat <- t(v) %*% W %*% mT2^2 %*% rR
    Kta <- t(Kat)
    Kaprec <- t(v) %*% D %*% mT2 %*% vI
    Kpreca <- t(Kaprec)
    
    Kpp <- t(rP) %*% W %*% mT2^2 %*% rP
    Kpt <- t(rP) %*% W %*% mT2^2 %*% rR
    Ktp <- t(Kpt)
    Kpprec <- t(rP) %*% D %*% mT2 %*% vI 
    Kprecp <- t(Kpprec)
    
    Ktt <- t(rR) %*% W %*% mT2^2 %*% rR
    Ktprec <- t(rR) %*% D %*% mT2 %*% vI
    Kprect <- t(Ktprec)
    
    Kprecprec <- sum(N)
    
    K <- -rbind(
      cbind(Koo,o,op,ot,o),
      cbind(t(o),Kaa,Kap,Kat,Kaprec),
      cbind(t(op),Kpa,Kpp,Kpt,Kpprec),
      cbind(t(ot),Kta,Ktp,Ktt,Ktprec),
      cbind(t(o),Kpreca,Kprecp,Kprect,Kprecprec)
    )

    z$K <- K

    z$K <- K
    z$vcov <- try(solve(z$K))
    stderror<-try(sqrt(diag(z$vcov)))
    
    z$stderror <- c(stderror[-(1:2)], stderror[1:2])

    if(anyNA(stderror))
    {
      z$Knumerica <- -gradient(escore, opt$par)
      z$vcov <- try(solve(z$Knumerica))
      z$stderror <- try(sqrt(diag(z$vcov)))
    }
    
    ynew_prev <- c(linkfun(ys),rep(NA,h1))
    errorhat_prev <- c(errorhat,rep(0,h1))
    y_prev <- c(z$fitted, rep(NA,h1))
    
    lambda_prev<- c(linkinv(etahat1), rep(NA,h1))
    s_prev<- c(rep(NA,m),z$s,rep(NA,h1))
    y_s <- c(y, rep(NA,h1))
    etahat1_prev <- c(z$etahat1, rep(NA,h1))

    for(i in 1:h1)
    {
      s_prev[n+i-1] <- meus(yt=y_s[1:(n+i-1)], b=b, n=n,link)
      
      etahat1_prev[n+i] <- omega1 + omega2*s_prev[n+i-1]
      lambda_prev[n+i] <- linkinv(etahat1_prev[n+i])
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat_prev[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      y_s[n+i] <- y_prev[n+i]
    }
  }
  
  ##################################################################################
  # ######### ARMA model sem regressores para [0,1) - IKARMA
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==T) && n0!=0 && n1==0)
  {
    reg <- c(rep(0,(p1+q1+1)), 10, 0.2, 0.1) 
    b <- 0
    reg <- c(linkfun(y_bar),rep(0,p1+q1), 0.2*prec, linkfun(delta0), 0)
    
    loglik <- function(d) 
    {
      alpha <- d[1]
      phi = d[2:(p1+1)] 
      theta = d[(p1+2):(p1+q1+1)]
      prec <- d[p1+q1+2] 
      omega1 <- d[p1+q1+3]

      error<-rep(0,n) 
      eta1<- eta2<- s <- rep(NA,n)

      for(i in (m+1):n)
      {
        s[i-1] <- meus(yt=y[1:(i-1)], b=b, n=n,link)
        eta1[i]<- omega1 + omega2*s[i-1]
        eta2[i] <- alpha + (phi%*%linkfun(ys[i-ar])) + (theta%*%error[i-ma])
        error[i] <- linkfun(ys[i])-eta2[i]
      }
      lambda <- linkinv(eta1[(m+1):n]) 
      mu <- linkinv(eta2[(m+1):n]) 
      
      l1 <- ifelse( (y1==0) | (y1==1), log(lambda),0) +
        ifelse( (y1==0) | (y1==1), 0,log(1-lambda)) 
      
      l2 <- ifelse((y1==0) |(y1==1), 0, 
                   log(prec) + 
                     ( log(log(0.5)/log(1-mu^prec)) ) +
                     (prec-1)*(log(y1)) +
                     ( (log(0.5)/log(1-mu^prec))-1)* ( log(1-y1^prec) ) )
      
      l3a <- ifelse( (y1==1),log(b),0)
      l3b <- ifelse( (y1==0),log(1-b),0)
      l3 <- l3a + l3b

      sum(l1 + l2+ l3)
    } 
    
    escore <- function(w)
    {
      alpha <- w[1]
      phi = w[2:(p1+1)]
      theta = w[(p1+2):(p1+q1+1)]
      prec <- w[p1+q1+2] 
      omega1 <- w[p1+q1+3]
      omega2 <- w[p1+q1+4]

      error<-rep(0,n) 
      eta1<- eta2<- s <- rep(NA,n)

      for(i in (m+1):n)
      {
        s[i-1] <- meus(yt=y[1:(i-1)], b=b, n=n,link)
        eta1[i]<- omega1 + omega2*s[i-1]
        eta2[i] <- alpha + (phi%*%linkfun(ys[i-ar])) + (theta%*%error[i-ma])
        error[i] <- linkfun(ys[i])-eta2[i]
      }
      
      lambda <- linkinv(eta1[(m+1):n])
      mu <- linkinv(eta2[(m+1):n])
      y1<-y[(m+1):n]
      s <- s[m:(n-1)]
      
      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        R[i,] <- error[i+m-ma]
      }

      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      }
      
      a <- ifelse( (y1 == 0) , (1)/(lambda), (-1)/((1-lambda)))
      a <- ifelse( (y1 == 1) , (1)/(lambda), a)
      
      I <- as.matrix(rep(1,n-m))
      
      delta <- log(0.5)/log(1-mu^prec)
      
      ci <- ifelse((y1 == 0) | (y1 == 1), 0, 
                   ( (mu^(prec-1)) / ((1-mu^prec)*log(1-mu^prec)) )*( delta*log(1-y1^prec)+1 ))
      
      c <- prec*ci

      nu <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rR <- deta.dtheta[(m+1):n,]
      
      mT1 <- diag(mu.eta(eta1[(m+1):n]))
      mT2 <- diag(mu.eta(eta2[(m+1):n]))
      
      phi1 <- 1/prec
      phi2 <- log(y1)
      phi3 <- ci*mu*log(mu)
      phi4 <- delta-1
      phi5 <- ( (y1^prec)*log(y1) ) / (1-y1^prec)
      dphi <- as.vector(phi1+phi2+phi3-(phi4*phi5))
      dphi1 <- ifelse( (y1==0)|(y1==1), 0,dphi) 
      
      Ualpha <- t(nu) %*% mT2 %*% c
      Uphi <-   t(rP) %*% mT2 %*% c
      Utheta <- t(rR) %*% mT2 %*% c
      Uprec <- sum(dphi1)
      Uomega1 <- t(a) %*% mT1 %*% I
      Uomega2 <- t(a) %*% mT1 %*% s

      rval <- c(Ualpha, Uphi, Utheta, Uprec, Uomega1, Uomega2)
    }
    names_par <- c("alpha", names_phi, names_theta, "precision","omega1", "omega2")
    
    
    opt <- optim(reg, loglik, escore, 
                 method = "BFGS", 
                 control = list(fnscale = -1, reltol = 1e-12,maxit=200))

    if (opt$conv != 0) 
      warning("FUNCTION DID NOT CONVERGE!")
    
    z$conv <- opt$conv
    coef <- (opt$par)
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    prec <- coef[p1+q1+2] 
    omega1 <- coef[p1+q1+3]
    omega2 <- coef[p1+q1+4]
    
    z$alpha<- alpha
    z$phi <- phi
    z$theta <- theta
    z$prec <- prec
    z$omega1 <- omega1 
    z$omega2 <- omega2 
    
    errorhat<-rep(0,n) 
    etahat1 <- etahat2 <- s <- rep(NA,n)
    
    for(i in (m+1):n)
    {
      s[i-1] <- meus(yt=y[1:(i-1)], b=b, n=n,link)
      etahat1[i] <- omega1 + omega2*s[i-1]
      etahat2[i] <-alpha + (phi%*%linkfun(ys[i-ar])) + (theta%*%errorhat[i-ma])
      errorhat[i] <- linkfun(ys[i])-etahat2[i]
    }
    
    lambdahat <- linkinv(etahat1[(m+1):n])
    muhat <- linkinv(etahat2[(m+1):n])
    y1<-y[(m+1):n]
    s <- s[m:(n-1)]
    #print(s)
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat1 <- etahat1
    z$etahat2 <- etahat2
    z$errorhat <- errorhat
    z$s <- s
    
    z$mu <- muhat
    z$lambda <- lambdahat
    z$pzero <- lambdahat*(1-b)
    z$pum <- lambdahat*(b)
    
    Q <- cbind(1,s)
    
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      R[i,] <- errorhat[i+m-ma]
    }
    
    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    
    for(i in (m+1):n)
    {
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
    }
    
    L <- diag(-1/((lambdahat)*(1-lambdahat)))
    
    v <- matrix(as.vector(deta.dalpha[(m+1):n]),ncol=1)
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    
    mT1 <- diag(mu.eta(etahat1[(m+1):n]))
    mT2 <- diag(mu.eta(etahat2[(m+1):n]))
    
    delta <- log(0.5)/log(1-muhat^prec)
    
    ci <- ifelse((y1 == 0) | (y1 == 1), 0, 
                 ((muhat^(prec-1))/((1-muhat^prec)*log(1-muhat^prec)))*(delta*log(1-y1^prec)+1))
    Vc <- prec*ci
    
    zeta1 <- (muhat^(prec-2))/((1-muhat^prec) * (log(1-muhat^prec)))
    zeta2 <- (muhat^(2*prec-2))/((1-muhat^prec)^2 * (log(1-muhat^prec))^2)
    
    W <- diag(as.vector((lambdahat -1)*  prec^2 * zeta2))
    vI <- matrix(rep(1,(n-m)),ncol=1)
    euler <- 0.5772156649
    
    d <- -prec*muhat*log(muhat)*zeta2-prec*muhat*delta*zeta1*((1-digamma(delta+1)-euler)/((delta-1)*prec))
    D <- diag(as.vector(d))
    
    N1 <- 1/(prec^2) 
    N2<- delta*(muhat^2)*zeta2*(log(muhat)^2)*log(1-y1^prec)
    N3<- 2*delta*(muhat^2)*log(muhat)*zeta1*(( (y1^prec) * log(y1))/(1-y1^prec))
    N4<- (delta-1)*( ((y1^prec)*log(y1)^2)/((1-y1^prec)^2) )   
    N5<- ci*muhat*(log(muhat)^2)*( (zeta1*(muhat^2)) +(1/(1-muhat^prec)) )
    
    N<- ifelse((y1 == 0), 0,-N1+N2-N3-N4+N5)
    N<- ifelse((y1 == 1), 0, N)
    
    o = t(matrix(data = 0, nrow = 1, ncol = 2))
    op = t(matrix(data = 0, nrow = p1, ncol = 2))
    ot = t(matrix(data = 0, nrow = q1, ncol = 2))
    
    Koo <- t(Q) %*% L %*% mT1^2 %*% Q
    
    Kaa <- t(v) %*% W %*% mT2^2 %*% v
    Kap <- t(v) %*% W %*% mT2^2 %*% rP
    Kpa <- t(Kap)
    Kat <- t(v) %*% W %*% mT2^2 %*% rR
    Kta <- t(Kat)
    Kaprec <- t(v) %*% D %*% mT2 %*% vI
    Kpreca <- t(Kaprec)
    
    Kpp <- t(rP) %*% W %*% mT2^2 %*% rP
    Kpt <- t(rP) %*% W %*% mT2^2 %*% rR
    Ktp <- t(Kpt)
    Kpprec <- t(rP) %*% D %*% mT2 %*% vI 
    Kprecp <- t(Kpprec)
    
    Ktt <- t(rR) %*% W %*% mT2^2 %*% rR
    Ktprec <- t(rR) %*% D %*% mT2 %*% vI
    Kprect <- t(Ktprec)
    
    Kprecprec <- sum(N)
    
    K <- -rbind(
      cbind(Koo,o,op,ot,o),
      cbind(t(o),Kaa,Kap,Kat,Kaprec),
      cbind(t(op),Kpa,Kpp,Kpt,Kpprec),
      cbind(t(ot),Kta,Ktp,Ktt,Ktprec),
      cbind(t(o),Kpreca,Kprecp,Kprect,Kprecprec)
    )

    z$K <- K
    
    z$K <- K
    z$vcov <- try(solve(z$K))
    stderror<-try(sqrt(diag(z$vcov)))
    
    z$stderror <- c(stderror[-(1:2)], stderror[1:2])
    if(anyNA(stderror))
    {
      z$Knumerica <- -gradient(escore, opt$par)
      z$vcov <- try(solve(z$Knumerica))
      z$stderror <- try(sqrt(diag(z$vcov)))
    }

    ynew_prev <- c(linkfun(ys),rep(NA,h1))
    errorhat_prev <- c(errorhat,rep(0,h1))
    y_prev <- c(z$fitted, rep(NA,h1))
    
    lambda_prev<- c(linkinv(etahat1), rep(NA,h1))
    s_prev<- c(rep(NA,m),z$s,rep(NA,h1))
    y_s <- c(y, rep(NA,h1))
    etahat1_prev <- c(z$etahat1, rep(NA,h1))
    
    for(i in 1:h1)
    {
      ##
      s_prev[n+i-1] <- meus(yt=y_s[1:(n+i-1)], b=b, n=n,link)
      
      etahat1_prev[n+i] <- omega1 + omega2*s_prev[n+i-1]
      lambda_prev[n+i] <- linkinv(etahat1_prev[n+i])
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat_prev[n+i-ma])
      
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      
      y_s[n+i] <- y_prev[n+i]
    }
  }
    
    
  
  ##################################################################################
  # ######### IKARMAX model com regressores para (0,1]
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==F) && n0==0 && n1!=0)
  { 
    beta1<- mqo[(p1+2):length(mqo)]
    b <- 1
    reg <- c(linkfun(y_bar),rep(0,p1+q1), 0.2*prec, linkfun(delta1), 0, rep(0,length(beta1)))
    
    loglik <- function(d) 
    {
      alpha <- d[1]
      phi = d[2:(p1+1)] 
      theta = d[(p1+2):(p1+q1+1)]
      prec <- d[p1+q1+2] 
      omega1 <- d[p1+q1+3]
      omega2 <- d[p1+q1+4]
      beta <- d[(p1+q1+5):length(d)]
      
      error<-rep(0,n) 
      eta1<- eta2<- s <- rep(NA,n)
      
      for(i in (m+1):n)
      {
        s[i-1] <- meus(yt=y[1:(i-1)], b=b, n=n,link)
        eta1[i]<- omega1 + omega2*s[i-1]
        eta2[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(linkfun(ys[i-ar])- X[i-ar,]%*%as.matrix(beta))) + (theta%*%error[i-ma])
        error[i] <- linkfun(ys[i])-eta2[i]
        
      }
      lambda <- linkinv(eta1[(m+1):n]) 
      mu <- linkinv(eta2[(m+1):n]) 
      
      l1 <- ifelse( (y1==0) | (y1==1), log(lambda),0) +
        ifelse( (y1==0) | (y1==1), 0,log(1-lambda))
      
      l2 <- ifelse((y1==0) |(y1==1), 0, 
                   log(prec) + 
                     ( log(log(0.5)/log(1-mu^prec)) ) +
                     (prec-1)*(log(y1)) +
                     ( (log(0.5)/log(1-mu^prec))-1)* ( log(1-y1^prec) ) )

      l3a <- ifelse( (y1==1),log(b),0)
      l3b <- ifelse( (y1==0),log(1-b),0)
      l3 <- l3a + l3b
      
      sum(l1 + l2+ l3)
      
    } 
    
    escore <- function(w)
    {
      alpha <- w[1]
      phi = w[2:(p1+1)]
      theta = w[(p1+2):(p1+q1+1)]
      prec <- w[p1+q1+2] 
      omega1 <- w[p1+q1+3]
      omega2 <- w[p1+q1+4]
      beta <- w[(p1+q1+5):length(w)]
      
      error<-rep(0,n) 
      eta1<- eta2<- s <- rep(NA,n)
      
      for(i in (m+1):n)
      {
        s[i-1] <- meus(yt=y[1:(i-1)], b=b, n=n,link)
        eta1[i]<- omega1 + omega2*s[i-1]
        eta2[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(linkfun(ys[i-ar])- X[i-ar,]%*%as.matrix(beta))) + (theta%*%error[i-ma])
        error[i] <- linkfun(ys[i])-eta2[i]
        
      }
      lambda <- linkinv(eta1[(m+1):n])
      mu <- linkinv(eta2[(m+1):n])
      y1<-y[(m+1):n]
      s <- s[m:(n-1)]
      
      P <- matrix(0,nrow=(n-m),ncol=p1)
      
      for(i in 1:(n-m))
      {
        P[i,] <- linkfun(ys[i+m-ar])- X[i+m-ar,]%*%as.matrix(beta)
      }
      
      R <- matrix(0,nrow=(n-m),ncol=q1)
      
      for(i in 1:(n-m))
      {
        R[i,]<-error[i+m-ma]
      }

      k1<- length(beta)
      
      M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
      for(i in 1:(n-m))
      {
        for(j in 1:k1)
          M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])
      }
      
      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dbeta<- matrix(0, ncol=k1,nrow=n)
 
      for(i in (m+1):n)
      {
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        deta.dbeta[i,]<- M[(i-m),] - theta%*% deta.dbeta[i-ma,]
      }
      
      nu <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rR <- deta.dtheta[(m+1):n,]
      rM <- deta.dbeta[(m+1):n,]

      a <- ifelse( (y1 == 0) , (1)/(lambda), (-1)/((1-lambda)))
      a <- ifelse( (y1 == 1) , (1)/(lambda), a)
      
      I <- as.matrix(rep(1,n-m))
      
      delta <- log(0.5)/log(1-mu^prec)
      
      ci <- ifelse((y1 == 0) | (y1 == 1), 0, 
                   ( (mu^(prec-1)) / ((1-mu^prec)*log(1-mu^prec)) )*( delta*log(1-y1^prec)+1 ))
      c <- prec*ci
      
      
      mT1 <- diag(mu.eta(eta1[(m+1):n]))
      mT2 <- diag(mu.eta(eta2[(m+1):n]))
      
      phi1 <- 1/prec
      phi2 <- log(y1)
      phi3 <- ci*mu*log(mu)
      phi4 <- delta-1
      phi5 <- ((y1^prec)*log(y1))/(1-y1^prec)
      dphi <- as.vector(phi1+phi2+phi3-(phi4*phi5)) 
      dphi1 <- ifelse( (y1==0)|(y1==1), 0,dphi) 
      
      Ualpha <- t(nu) %*% mT2 %*% c
      Uphi <-   t(rP) %*% mT2 %*% c
      Utheta <- t(rR) %*% mT2 %*% c
      Uprec <- sum(dphi1)
      Uomega1 <- t(a) %*% mT1 %*% I
      Uomega2 <- t(a) %*% mT1 %*% s
      Ubeta <- t(rM) %*% mT2 %*% c
      
      
      rval <- c(Ualpha, Uphi, Utheta, Uprec, Uomega1, Uomega2, Ubeta)
      
    }
    
    names_par <- c("alpha",names_phi,names_theta,"precision","omega1","omega2",names_beta)
    
    opt <- optim(reg, loglik, escore, 
                 method = "BFGS", 
                 control = list(fnscale = -1))
    
    if (opt$conv != 0) 
      warning("FUNCTION DID NOT CONVERGE!")
    
    z$conv <- opt$conv
    coef <- (opt$par)
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    prec <- coef[p1+q1+2] 
    omega1 <- coef[p1+q1+3]
    omega2 <- coef[p1+q1+4]
    beta <- coef[(p1+q1+5):length(coef)]
    
    z$alpha<- alpha
    z$phi <- phi
    z$theta <- theta
    z$prec <- prec
    z$omega1 <- omega1 
    z$omega2 <- omega2 
    z$beta <- beta

    errorhat<-rep(0,n) 
    etahat1 <- etahat2 <- s <- rep(NA,n)
    
    
    for(i in (m+1):n)
    {
      s[i-1] <- meus(yt=y[1:(i-1)], b=b, n=n,link)
      etahat1[i] <- omega1 + omega2*s[i-1]
      etahat2[i] <-alpha + X[i,]%*%as.matrix(beta) + (phi%*%(linkfun(ys[i-ar]) - X[i-ar,]%*%as.matrix(beta))) + (theta%*%errorhat[i-ma])
      errorhat[i] <- linkfun(ys[i])-etahat2[i]
      
    }
    lambdahat <- linkinv(etahat1[(m+1):n])
    muhat <- linkinv(etahat2[(m+1):n])
    y1<-y[(m+1):n]
    s <- s[m:(n-1)]
    
    z$mu <- muhat
    z$lambda <- lambdahat
    z$pzero <- lambdahat*(1-b)
    z$pum <- lambdahat*(b)
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat1 <- etahat1
    z$etahat2 <- etahat2
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    z$s <- s

    Q <- cbind(1,s)
    
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      R[i,] <- errorhat[i+m-ma]
    }
    
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    for(i in 1:(n-m))
    {
      P[i,] <- linkfun(ys[i+m-ar]) - X[i+m-ar,]%*%as.matrix(beta)
    }
    
    k1<- length(beta)
    
    M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
    for(i in 1:(n-m))
    {
      for(j in 1:k1)
        M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])
    }

    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dbeta<- matrix(0, ncol=k1,nrow=n)
    
    for(i in (m+1):n)
    {
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
    }
    
    L <- diag(-1/((lambdahat)*(1-lambdahat)))
    
    v <- matrix(as.vector(deta.dalpha[(m+1):n]),ncol=1)
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    rM <- deta.dbeta[(m+1):n,]
    
    mT1 <- diag(mu.eta(etahat1[(m+1):n]))
    mT2 <- diag(mu.eta(etahat2[(m+1):n]))
    
    delta <- log(0.5)/log(1-muhat^prec)
    
    ci <- ifelse((y1 == 0) | (y1 == 1), 0, 
                 ((muhat^(prec-1))/((1-muhat^prec)*log(1-muhat^prec)))*(delta*log(1-y1^prec)+1) )
    Vc <- prec*ci
    
    zeta1 <- (muhat^(prec-2))/((1-muhat^prec) * (log(1-muhat^prec)))
    zeta2 <- (muhat^(2*prec-2))/((1-muhat^prec)^2 * (log(1-muhat^prec))^2)
    
    W <- diag(as.vector((lambdahat -1)*  prec^2 * zeta2))
    vI <- matrix(rep(1,(n-m)),ncol=1)
    
    euler <- 0.5772156649
    
    d <- -prec*muhat*log(muhat)*zeta2-prec*muhat*delta*zeta1*((1-digamma(delta+1)-euler)/((delta-1)*prec))
    D <- diag(as.vector(d))
    N1 <- 1/(prec^2) 
    N2<- delta*(muhat^2)*zeta2*(log(muhat)^2)*log(1-y1^prec)
    N3<- 2*delta*(muhat^2)*log(muhat)*zeta1*(( (y1^prec) * log(y1))/(1-y1^prec))
    N4<- (delta-1)*((y1^prec)*log(y1)^2)/((1-y1^prec)^2)   
    N5<- Vc*muhat*(log(muhat)^2)*(zeta2*(muhat^2)+1/(1-muhat^prec))
    
    N<- ifelse((y1 == 0), 0,-N1+N2-N3-N4+N5)
    N<- ifelse((y1 == 1), 0, N)
    
    Koo <- t(Q) %*% L %*% mT1^2 %*% Q
    
    Kaa <- t(v) %*% W %*% mT2^2 %*% v
    Kap <- t(v) %*% W %*% mT2^2 %*% rP
    Kpa <- t(Kap)
    Kat <- t(v) %*% W %*% mT2^2 %*% rR
    Kta <- t(Kat)
    Kaprec <- t(v) %*% D %*% mT2 %*% vI
    Kpreca <- t(Kaprec)
    
    Kpp <- t(rP) %*% W %*% mT2^2 %*% rP
    Kpt <- t(rP) %*% W %*% mT2^2 %*% rR
    Ktp <- t(Kpt)
    Kpprec <- t(rP) %*% D %*% mT2 %*% vI 
    Kprecp <- t(Kpprec)
    
    Ktt <- t(rR) %*% W %*% mT2^2 %*% rR
    Ktprec <- t(rR) %*% D %*% mT2 %*% vI
    Kprect <- t(Ktprec)
    
    Kprecprec <- sum(N)
    
    Kab <- t(v) %*% W %*% mT2^2 %*% rM
    Kba <- t(Kab)
    Kbb <- t(rM) %*% W %*% mT2^2 %*% rM
    Kpb <- t(rP) %*% W %*% mT2^2 %*% rM
    Kbp <- t(Kpb)
    Ktb <- t(rR) %*% W %*% mT2^2 %*% rM
    Kbt <- t(Ktb)
    Kbprec <- t(rM) %*% D %*% mT2 %*% vI
    Kprecb <- t(Kbprec)
    
    o = t(matrix(data = 0, nrow = 1, ncol = 2))
    op = t(matrix(data = 0, nrow = p1, ncol = 2))
    ot = t(matrix(data = 0, nrow = q1, ncol = 2))
    ob = t(matrix(data = 0, nrow = k1, ncol = 2))
    
    K <- -rbind(
      cbind(Koo,o,op,ot,o,ob),
      cbind(t(o),Kaa,Kap,Kat,Kaprec,Kab),
      cbind(t(op),Kpa,Kpp,Kpt,Kpprec,Kpb),
      cbind(t(ot),Kta,Ktp,Ktt,Ktprec,Ktb),
      cbind(t(o),Kpreca,Kprecp,Kprect,Kprecprec,Kprecb),
      cbind(t(ob),Kba,Kbp,Kbt,Kbprec,Kbb)
    )
    
    z$K <- K
    z$vcov <- try(solve(z$K))
    stderror<-try(sqrt(diag(z$vcov)))
    
    z$stderror <- c(stderror[-c(1,2,(p1+q1+5):ncol(K))], stderror[c(1,2,(p1+q1+5):ncol(K))])
    if(anyNA(stderror))
    {
      z$Knumerica <- -gradient(escore, opt$par)
      z$vcov <- try(solve(z$Knumerica))
      z$stderror <- try(sqrt(diag(z$vcov)))
    }

    ynew_prev <- c(linkfun(ys),rep(NA,h1))
    errorhat_prev <- c(errorhat,rep(0,h1))
    y_prev <- c(z$fitted, rep(NA,h1))
    X_prev<- rbind(X,X_hat)

    lambda_prev<- c(linkinv(etahat1), rep(NA,h1))
    s_prev<- c(rep(NA,m),z$s,rep(NA,h1))
    y_s <- c(y, rep(NA,h1))
    etahat1_prev <- c(z$etahat1, rep(NA,h1))

    for(i in 1:h1)
    {
      s_prev[n+i-1] <- meus(yt=y_s[1:(n+i-1)], b=b, n=n,link)

      etahat1_prev[n+i] <- omega1 + omega2*s_prev[n+i-1]
      lambda_prev[n+i] <- linkinv(etahat1_prev[n+i])
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar] -X_prev[n+i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat_prev[n+i-ma])

      y_prev[n+i] <- linkinv(ynew_prev[n+i])

      y_s[n+i] <- y_prev[n+i]
    }
  
  } 
  
  
  ##################################################################################
  # ######### IKARMAX model com regressores para [0,1)
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==F) && n0!=0 && n1==0)
  { 
    
    b <- 0
    beta1<- mqo[(p1+2):length(mqo)]
    reg <- c(linkfun(y_bar),rep(0,p1+q1), 0.2*prec, linkfun(delta0), 0, rep(0,length(beta1)) )

    loglik <- function(d) 
    {
      alpha <- d[1]
      phi <- d[2:(p1+1)] 
      theta <- d[(p1+2):(p1+q1+1)]
      prec <- d[p1+q1+2] 
      omega1 <- d[p1+q1+3]
      omega2 <- d[p1+q1+4]
      beta <- d[(p1+q1+5):length(d)]
      
      error<-rep(0,n) #
      eta1<- eta2<- s <- rep(NA,n)
      
      for(i in (m+1):n)
      {
        s[i-1] <- meus(yt=y[1:(i-1)], b=b, n=n,link)
        eta1[i]<- omega1 + omega2*s[i-1]
        eta2[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(linkfun(ys[i-ar])- X[i-ar,]%*%as.matrix(beta))) + (theta%*%error[i-ma])
        error[i] <- linkfun(ys[i])-eta2[i]
        
      }
      lambda <- linkinv(eta1[(m+1):n]) 
      mu <- linkinv(eta2[(m+1):n])
      
      
      l1 <- ifelse( (y1==0) | (y1==1), log(lambda),0) +
            ifelse( (y1==0) | (y1==1), 0,log(1-lambda))
      
      l2 <- ifelse((y1==0) |(y1==1), 0, 
                   log(prec) + 
                     ( log(log(0.5)/log(1-mu^prec)) ) +
                     (prec-1)*(log(y1)) +
                     ( (log(0.5)/log(1-mu^prec))-1)* ( log(1-y1^prec) ) )
      

      l3a <- ifelse( (y1==1),log(b),0)
      l3b <- ifelse( (y1==0),log(1-b),0)
      l3 <- l3a + l3b
      
      sum(l1 + l2+ l3)
      
    } 
    
    escore <- function(w)
    {
      alpha <- w[1]
      phi = w[2:(p1+1)]
      theta = w[(p1+2):(p1+q1+1)]
      prec <- w[p1+q1+2] 
      omega1 <- w[p1+q1+3]
      omega2 <- w[p1+q1+4]
      beta <- w[(p1+q1+5):length(w)]
      
      error<-rep(0,n) 
      eta1<- eta2<- s <- rep(NA,n)
      
      
      for(i in (m+1):n)
      {
        s[i-1] <- meus(yt=y[1:(i-1)], b=b, n=n,link)
        eta1[i]<- omega1 + omega2*s[i-1]
        eta2[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(linkfun(ys[i-ar])- X[i-ar,]%*%as.matrix(beta))) + (theta%*%error[i-ma])
        error[i] <- linkfun(ys[i])-eta2[i]
        
      }
      lambda <- linkinv(eta1[(m+1):n])
      mu <- linkinv(eta2[(m+1):n])
      y1<-y[(m+1):n]
      s <- s[m:(n-1)]
      
      P <- matrix(0,nrow=(n-m),ncol=p1)
      
      for(i in 1:(n-m))
      {
        P[i,] <- linkfun(ys[i+m-ar])- X[i+m-ar,]%*%as.matrix(beta)
      }
      
      R <- matrix(0,nrow=(n-m),ncol=q1)
      
      for(i in 1:(n-m))
      {
        R[i,]<-error[i+m-ma]
      }
      
      k1<- length(beta)
      
      M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
      for(i in 1:(n-m))
      {
        for(j in 1:k1)
          M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])
      }
      
      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dbeta<- matrix(0, ncol=k1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        deta.dbeta[i,]<- M[(i-m),] - theta%*% deta.dbeta[i-ma,]
      }
      
      nu <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rR <- deta.dtheta[(m+1):n,]
      rM <- deta.dbeta[(m+1):n,]
     
      a <- ifelse( (y1 == 0) , (1)/(lambda), (-1)/((1-lambda)))
      a <- ifelse( (y1 == 1) , (1)/(lambda), a)
      
      I <- as.matrix(rep(1,n-m))
      
      delta <- log(0.5)/log(1-mu^prec)
      
      ci <- ifelse((y1 == 0) | (y1 == 1), 0, 
                   ( (mu^(prec-1)) / ((1-mu^prec)*log(1-mu^prec)) )*( delta*log(1-y1^prec)+1 ))
      
      c <- prec*ci
      
      
      mT1 <- diag(mu.eta(eta1[(m+1):n]))
      mT2 <- diag(mu.eta(eta2[(m+1):n]))
      
      phi1 <- 1/prec
      phi2 <- log(y1)
      phi3 <- ci*mu*log(mu)
      phi4 <- delta-1
      phi5 <- ((y1^prec)*log(y1))/(1-y1^prec)
      dphi <- as.vector(phi1+phi2+phi3-(phi4*phi5)) 
      dphi1 <- ifelse( (y1==0)|(y1==1), 0,dphi) 
      
      Ualpha <- t(nu) %*% mT2 %*% c
      Uphi <-   t(rP) %*% mT2 %*% c
      Utheta <- t(rR) %*% mT2 %*% c
      Uprec <- sum(dphi1)
      Uomega1 <- t(a) %*% mT1 %*% I
      Uomega2 <- t(a) %*% mT1 %*% s
      Ubeta <- t(rM) %*% mT2 %*% c
      
      
      rval <- c(Ualpha, Uphi, Utheta, Uprec, Uomega1, Uomega2, Ubeta)
      
    }
    
    names_par <- c("alpha",names_phi,names_theta,"precision","omega1","omega2",names_beta)
    
    opt <- optim(reg, loglik, escore, 
                 method = "BFGS", 
                 control = list(fnscale = -1))
    
    if (opt$conv != 0) 
    warning("FUNCTION DID NOT CONVERGE!")
    
    
    z$conv <- opt$conv
    coef <- (opt$par)
    names(coef)<-names_par
    z$coeff <- coef
    
   
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    prec <- coef[p1+q1+2] 
    omega1 <- coef[p1+q1+3]
    omega2 <- coef[p1+q1+4]
    beta <- coef[(p1+q1+5):length(coef)]
    
    z$alpha<- alpha
    z$phi <- phi
    z$theta <- theta
    z$prec <- prec
    z$omega1 <- omega1 
    z$omega2 <- omega2 
    z$beta <- beta
    #z$b <- b
    
    errorhat<-rep(0,n) 
    etahat1 <- etahat2 <- s <- rep(NA,n)
    
    
    for(i in (m+1):n)
    {
      s[i-1] <- meus(yt=y[1:(i-1)], b=b, n=n,link)
      etahat1[i] <- omega1 + omega2*s[i-1]
      etahat2[i] <-alpha + X[i,]%*%as.matrix(beta) + (phi%*%(linkfun(ys[i-ar]) - X[i-ar,]%*%as.matrix(beta))) + (theta%*%errorhat[i-ma])
      errorhat[i] <- linkfun(ys[i])-etahat2[i]
      
    }
    lambdahat <- linkinv(etahat1[(m+1):n])
    muhat <- linkinv(etahat2[(m+1):n])
    y1<-y[(m+1):n]
    s <- s[m:(n-1)]
    
    z$mu <- muhat
    z$lambda <- lambdahat
    z$pzero <- lambdahat*(1-b)
    z$pum <- lambdahat*(b)
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat1 <- etahat1
    z$etahat2 <- etahat2
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    z$s <- s
    
    Q <- cbind(1,s)
    
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      R[i,] <- errorhat[i+m-ma]
    }
    
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    for(i in 1:(n-m))
    {
      P[i,] <- linkfun(ys[i+m-ar]) - X[i+m-ar,]%*%as.matrix(beta)
    }
    
    k1<- length(beta)

    M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
    for(i in 1:(n-m))
    {
      for(j in 1:k1)
        M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])
    }
    
    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dbeta<- matrix(0, ncol=k1,nrow=n)
    
    for(i in (m+1):n)
    {
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
    }
    
    L <- diag(-1/((lambdahat)*(1-lambdahat)))
    
    v <- matrix(as.vector(deta.dalpha[(m+1):n]),ncol=1)
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    rM <- deta.dbeta[(m+1):n,]
    
    mT1 <- diag(mu.eta(etahat1[(m+1):n]))
    mT2 <- diag(mu.eta(etahat2[(m+1):n]))
    
    delta <- log(0.5)/log(1-muhat^prec)
    
    ci <- ifelse((y1 == 0) | (y1 == 1), 0, 
                 ((muhat^(prec-1))/((1-muhat^prec)*log(1-muhat^prec)))*(delta*log(1-y1^prec)+1) )
    
    Vc <- prec*ci
    
    
    zeta1 <- (muhat^(prec-2))/((1-muhat^prec) * (log(1-muhat^prec)))
    zeta2 <- (muhat^(2*prec-2))/((1-muhat^prec)^2 * (log(1-muhat^prec))^2)
    
    W <- diag(as.vector((lambdahat -1)*  prec^2 * zeta2))
    vI <- matrix(rep(1,(n-m)),ncol=1)
    
    euler <- 0.5772156649
    
    d <- -prec*muhat*log(muhat)*zeta2-prec*muhat*delta*zeta1*((1-digamma(delta+1)-euler)/((delta-1)*prec))
    D <- diag(as.vector(d))
    N1 <- 1/(prec^2) 
    N2<- delta*(muhat^2)*zeta2*(log(muhat)^2)*log(1-y1^prec)
    N3<- 2*delta*(muhat^2)*log(muhat)*zeta1*(( (y1^prec) * log(y1))/(1-y1^prec))
    N4<- (delta-1)*((y1^prec)*log(y1)^2)/((1-y1^prec)^2)   
    N5<- Vc*muhat*(log(muhat)^2)*(zeta2*(muhat^2)+1/(1-muhat^prec))
    
    N<- ifelse((y1 == 0), 0,-N1+N2-N3-N4+N5)
    N<- ifelse((y1 == 1), 0, N)
    
    Koo <- t(Q) %*% L %*% mT1^2 %*% Q
    
    Kaa <- t(v) %*% W %*% mT2^2 %*% v
    Kap <- t(v) %*% W %*% mT2^2 %*% rP
    Kpa <- t(Kap)
    Kat <- t(v) %*% W %*% mT2^2 %*% rR
    Kta <- t(Kat)
    Kaprec <- t(v) %*% D %*% mT2 %*% vI
    Kpreca <- t(Kaprec)
    
    Kpp <- t(rP) %*% W %*% mT2^2 %*% rP
    Kpt <- t(rP) %*% W %*% mT2^2 %*% rR
    Ktp <- t(Kpt)
    Kpprec <- t(rP) %*% D %*% mT2 %*% vI 
    Kprecp <- t(Kpprec)
    
    Ktt <- t(rR) %*% W %*% mT2^2 %*% rR
    Ktprec <- t(rR) %*% D %*% mT2 %*% vI
    Kprect <- t(Ktprec)
    
    Kprecprec <- sum(N)
    
    Kab <- t(v) %*% W %*% mT2^2 %*% rM
    Kba <- t(Kab)
    Kbb <- t(rM) %*% W %*% mT2^2 %*% rM
    Kpb <- t(rP) %*% W %*% mT2^2 %*% rM
    Kbp <- t(Kpb)
    Ktb <- t(rR) %*% W %*% mT2^2 %*% rM
    Kbt <- t(Ktb)
    Kbprec <- t(rM) %*% D %*% mT2 %*% vI
    Kprecb <- t(Kbprec)
    
    o = t(matrix(data = 0, nrow = 1, ncol = 2))
    op = t(matrix(data = 0, nrow = p1, ncol = 2))
    ot = t(matrix(data = 0, nrow = q1, ncol = 2))
    ob = t(matrix(data = 0, nrow = k1, ncol = 2))
    
    
    K <- -rbind(
      cbind(Koo,o,op,ot,o,ob),
      cbind(t(o),Kaa,Kap,Kat,Kaprec,Kab),
      cbind(t(op),Kpa,Kpp,Kpt,Kpprec,Kpb),
      cbind(t(ot),Kta,Ktp,Ktt,Ktprec,Ktb),
      cbind(t(o),Kpreca,Kprecp,Kprect,Kprecprec,Kprecb),
      cbind(t(ob),Kba,Kbp,Kbt,Kbprec,Kbb)
    )

    z$K <- K
    z$vcov <- try(solve(z$K))
    stderror<-try(sqrt(diag(z$vcov)))
    
    z$stderror <- c(stderror[-c(1,2,(p1+q1+5):ncol(K))], stderror[c(1,2,(p1+q1+5):ncol(K))])
    if(anyNA(stderror))
    {
      z$Knumerica <- -gradient(escore, opt$par)
      z$vcov <- try(solve(z$Knumerica))
      z$stderror <- try(sqrt(diag(z$vcov)))
    }
    
    ynew_prev <- c(linkfun(ys),rep(NA,h1))
    errorhat_prev <- c(errorhat,rep(0,h1))
    y_prev <- c(z$fitted, rep(NA,h1))
    X_prev<- rbind(X,X_hat)
    
    lambda_prev<- c(linkinv(etahat1), rep(NA,h1))
    s_prev<- c(rep(NA,m),z$s,rep(NA,h1))
    y_s <- c(y, rep(NA,h1))
    etahat1_prev <- c(z$etahat1, rep(NA,h1))
    
    for(i in 1:h1)
    {
      s_prev[n+i-1] <- meus(yt=y_s[1:(n+i-1)], b=b, n=n,link)
      etahat1_prev[n+i] <- omega1 + omega2*s_prev[n+i-1]
      lambda_prev[n+i] <- linkinv(etahat1_prev[n+i])
      
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar] -X_prev[n+i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat_prev[n+i-ma])
      
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      y_s[n+i] <- y_prev[n+i]
      
    }
    
  }  
  
  
  z$serie <- y
  z$ikarma <- names_par
  z$forecast <- y_prev[(n+1):(n+h1)]
  z$y_prev <- ts(y_prev,start=start(y),frequency=frequency(y))
  
  
  z$fitted2<-ts(c(rep(NA,m),qikarma(rep(0.5,length(z$lambda)),z$lambda, b=b, z$mu,z$prec)),start=start(y),frequency=frequency(y)) 
  z$forecast2<-qikarma(rep(0.5,h1),lambda_prev[(n+1):(n+h1)], b=b, z$forecast,z$prec)
  z$y_prev2 <- ts(c(z$fitted2,z$forecast2),start=start(y),frequency=frequency(y)) 
    
  
  #######################################################################################
  #######################################################################################
  #######################################################################################

  # residuals
  res1 <- y-z$fitted
  vary <- ( (log(0.5)/log(1-z$fitted^prec)) * beta(1+2/prec, log(0.5)/log(1-z$fitted^prec))
                      - (log(0.5)/log(1-z$fitted^prec) * beta(1+1/prec, log(0.5)/log(1-z$fitted^prec)))^2 )
 
  z$resid1 <- (res1/sqrt(vary))[(m+1):n]
  
  l_tilde <- log(dikarma(y, lambdahat,b, y, z$prec))
  l_hat <- log(dikarma(y, lambdahat, b, z$fitted,z$prec))
 
  dt <- (l_tilde-l_hat)[(m+1):n]
  dt[which(dt<0)]<-0

  z$l_hat <- l_hat

  z$resid2 <- sign(y[(m+1):n]-z$fitted[(m+1):n])*sqrt(2*(dt)) 

  z$resid3 <- as.vector(qnorm(pikarma(y[(m+1):n], lambdahat, b, z$fitted[(m+1):n],z$prec)))

  # randomized quantile  residuals
  ai0 <- 0 
  bi0 <- z$lambda*(1-b)
  
  ai1 <- 1 - z$lambda*b
  bi1 <- 1
  
  ru0 <- runif(n,ai0,bi0)
  ru1 <- runif(n,ai1,bi1)
  
  ui <- pikarma(y1,z$lambda,b,z$mu,z$prec)
  ui <- ifelse(y1==0,ru0,ui)
  ui <- ifelse(y1==1,ru1,ui)
  
  z$resid4 = qnorm(ui)
  
  if(resid==1) residc <- z$resid1
  if(resid==2) residc <- z$resid2
  if(resid==3) residc <- z$resid3
  if(resid==4) residc <- z$resid4

  ############################################

  z$zstat <- abs(z$coef/z$stderror)
  z$pvalues <- 2*(1 - pnorm(z$zstat) )

  z$loglik <- opt$value
  z$counts <- as.numeric(opt$counts[1])

  if(any(is.na(X)==F))
  {
    z$k<- (p1+q1+2+length(beta))
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }else{
    z$k<- (p1+q1+2)
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }

  model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")

  z$model <- model_presentation
  z$link <- link

  z$escore <- escore(opt$par)
 
  ###################################################
  ######### GRAPHICS ################################
  
  if(diag>0)
  {
    
    print(model_presentation)
    print(" ",quote=F)
    print(c("Log-likelihood:",round(z$loglik,4)),quote=F)
    print(c("Number of iterations in BFGS optim:",z$counts),quote=F)
    print(c("AIC:",round(z$aic,4)," BIC:",round(z$bic,4)),quote=F)
    
    print("Residuals:",quote=F)
    print(summary(residc))
    
    t<-seq(-5,n+6,by=1)
    
    par(mfrow=c(1,1))
    par(mar=c(2.8, 2.7, 1.2, 1)) 
    par(mgp=c(1.7, 0.45, 0))
    plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
    lines(t,rep(-3,n+12),lty=2,col=1)
    lines(t,rep(3,n+12),lty=2,col=1)
    lines(t,rep(-2,n+12),lty=3,col=1)
    lines(t,rep(2,n+12),lty=3,col=1)

    
    densidade<-density(residc)
    plot(densidade,ylab="Density",main=" ")
    lines(densidade$x,dnorm(densidade$x),lty=2)
    legend("topleft",c("Density","Standard normal"),
           pt.bg="white", lty=c(1,2), bty="n")
    
    acf(residc,ylab="ACF",xlab="lag") 
    
    pacf(residc,ylab="PACF",xlab="lag") 
    
    max_r<- max(residc,na.rm=T)
    min_r<- min(residc,na.rm=T)
    qqnorm(residc, pch = "+",
           xlim=c(0.95*min_r,max_r*1.05),
           ylim=c(0.95*min_r,max_r*1.05),
           main="",xlab="Normal quantile",ylab="Empirical quantile")
    lines(c(-10,10),c(-10,10),lty=2)
    
    
    fim<-end(y)[1]+end(y)[2]/12
    
    par(mfrow=c(1,1))
    par(mar=c(2.8, 2.7, 1, 1)) 
    par(mgp=c(1.7, 0.45, 0))
    plot(y,type="l",ylab="Time series",xlab="Time")
    lines(z$fitted,col="red")
    
    w1<-3
    h1<-3
    
    if(diag>1)
    {
      pdf(file = "resid_v_ind.pdf",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
        lines(t,rep(-3,n+12),lty=2,col=1)
        lines(t,rep(3,n+12),lty=2,col=1)
        lines(t,rep(-2,n+12),lty=3,col=1)
        lines(t,rep(2,n+12),lty=3,col=1)
      }
      dev.off()
      
      
      pdf(file = "resid_density.pdf",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(1.5, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        
        plot(densidade,ylab="Density",main=" ",xlab=" ",ylim=c(0,1.15*max(densidade$y)))
        lines(densidade$x,dnorm(densidade$x),lty=2)
        legend("topleft",c("Fitted dansity","Standard normal"),
               pt.bg="white", lty=c(1,2), bty="n")
      }
      dev.off()
      
      pdf(file = "resid_FAC.pdf",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        acf(residc,ylab="ACF",xlab="lag") 
      }
      dev.off()
      
      pdf(file = "resid_FACP.pdf",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        pacf(residc,ylab="PACF",xlab="lag") 
      }
      dev.off()
      
      pdf(file = "qq_plot.pdf",width = w1, height = h1,family = "Times")
      {  
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        qqnorm(residc, pch = "+",
               xlim=c(0.95*min_r,max_r*1.05),
               ylim=c(0.95*min_r,max_r*1.05),
               main="",xlab="Normal quantile",ylab="Empirical quantile")
        lines(c(-10,10),c(-10,10),lty=2)
      }
      dev.off()
      
      pdf(file = "adjusted.pdf",width = 6, height = 4,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        plot(y,type="l",ylab="Data",xlab="Time")
        lines(z$fitted,col="red")
      }
      dev.off()
      
      pdf(file = "forecast.pdf",width = 6, height = 4,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        plot(y_prev,type="l",col="red",lty=2, ylim=c(min(y),max(y)),ylab="Data",xlab="Time")
        abline(v=fim,lty=2)
        lines(y)
      }
      dev.off()
    }    
  }
 
  return(z)
}







