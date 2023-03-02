# Implemented by Fabio M Bayer and Camila Malu da Rosa
# e-mail: bayer@ufsm.br

#######################################################################
######################## Monte Carlo simulation #######################
#######################################################################

source("ikarma.R")
source("simu_ikarma.R")

#set.seed(12345)
R<-10 # Monte Carlo replications
vn<-c(100,300,500,1000) # sample sizes

# parameter values
alpha<- -1
phi<- -0.45
theta<- 0.3
prec<- 5
omega1<- -3
omega2<- 0.5
b<-0

pa<-c(alpha,phi,theta,prec,omega1,omega2)

for(n in vn)
{
  TC<-rep(0,length(pa))
  coef_results<-c()
  results<-c()
  i<-1
  nc<-0 
  problemas <- 0 
  pinfl<-c()
  
  while (i <= R)
  {
    y<-simu.ikarma(n,alpha=alpha,phi=phi,theta=theta, prec=prec,
                   omega1=omega1, omega2=omega2, b=b) 
    
    if( ((sum(y==1)+sum(y==0))!=0)) 
    {
      fit1<-try(ikarma(y,ar=c(1),ma=c(1),diag=0)) 
      
      if((class(fit1) != "try-error") )
      {
        if((all(fit1$escore)<200) )
        {
          if(fit1$conv==0)
          {
            coef_results<-rbind(coef_results,fit1$coef)
            
            results<-rbind(results,fit1$pvalues)
            
            LI<- fit1$coef - qnorm(0.975)* fit1$stderror
            LS<- fit1$coef + qnorm(0.975)* fit1$stderror
            
            TC <- TC + ( (pa<LI) + (pa>LS))
            
            pinfl<-c(pinfl,100*((sum(y==1)+sum(y==0)))/n)
            
            i<-i+1
          }else{ 
            nc<-nc+1
            print(c("Non convergence",i,nc))
          }
        }else{problemas<-problemas+1}
      }else{problemas<-problemas+1}   
    }
  }
  
  m<-colMeans((coef_results),na.rm=T)
  mediana<-apply(coef_results,2,median)
  sd<-apply(coef_results,2,sd)
  bias<-m-pa
  rb<- 100*bias/pa
  tc <- 1-TC/R 
  
  M<-rbind(pa,m,mediana,sd,bias,rb, tc)
  colnames(M)<-c("alpha", "phi", "theta", "prec", "omega1","omega2")
  row.names(M)<-c("Parameters","Mean","Median","EP","Bias","RB", "TC")
  
  print(" ", quote=F)
  print(c("Sample size =",n),quote=F) 
  print(round(M,3))
  print(c("n=",n, " Non converges=",nc, " Numerical erros=",problemas," Mean of inflation=",mean(pinfl)),quote=F)
  
}

