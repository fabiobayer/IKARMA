# Implemented by Fabio M Bayer and Camila Malu da Rosa
# e-mail: bayer@ufsm.br

source("ikarma.R")

data<-read.table("ur-max-mensal-caxias.txt") 
attach(data)
y1<-x/100 # 
N<-length(y1)
y<-y1[1:(N-12)]
yout<-y1[(N-11):N]
y<-ts(y,start=c(2002,1),frequency=12) 
n<-length(y)
n

plot(y) 
acf(y,main="")
pacf(y,main="")

plot(decompose(y))
monthplot(y)

summary(y)
sd(y) 
CV <- 100*sd(y)/mean(y) 
CV

# proportion of inflation
sum(y==1)
sum(y==1)/n
sum(y==0)
sum(y==0)/n

h1<-12 

# fit without covariates
ikarma1<-ikarma(y,ar=c(1,2),ma=c(1,2),diag=1,resid=4,h=12) 

# tendency
m<-5
t <- (1+m):(n+m) # in sample
t_hat <- (n+1+m):(n+h1+m) # out of sample

# deterministic seasonality
S<-sin(2*pi*t/12) # in sample
S_hat<-sin(2*pi*t_hat/12) # out of sample

mX<-cbind(S) # in sample
mX_hat<-cbind(S_hat) # out of sample

# fit with covariate
ikarma2<-ikarma(y,ar=c(1),ma=c(1),diag=2, X = mX,X_hat = mX_hat,resid=4,h=h1) 

plot(ikarma2$resid4,type="p",pch="+")
acf(ikarma2$resid4)
pacf(ikarma2$resid4)

