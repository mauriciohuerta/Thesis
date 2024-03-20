library(Rfast)
library(forecast)

### Simulcaciones

n <- 400
t <- 20

p11 <- runif(1, -1, 1)
p12 <- runif(1, -1, 1)
p21 <- runif(1, -1, 1)
p22 <- runif(1, -1, 1)


q11 <- runif(1, -1, 1)
q12 <- runif(1, -1, 1)
q21 <- runif(1, -1, 1)
q22 <- runif(1, -1, 1)

s11 <- runif(1, 10, 100)
s22 <- runif(1, 20, 200)
s12 <- 0

p1=matrix(c(p11,p12,p21,p22),2,2)
sig=matrix(c(s11,s12,s12,s22),2,2)
th1=matrix(c(q11,q12,q21,q22),2,2)
m1=VARMAsim(n+t,arlags=c(1),malags=c(0),phi=p1,theta=th1,sigma=sig)
zt=m1$series

#Estimacion sin tiempo

data <- zt[1:n,]
simV <- zt[101:110,]

mu1 <- 200
mu2 <- 400

vest <- VARMA(data, p=1, q=1, include.mean = F)
vest$Sigma; sig

data[,1] <- data[,1]+mu1
data[,2] <- data[,2]+mu2

nest <- mvnorm.mle(data)
nest$mu

pred <- VARMApred(vest, h = t)

par(mfrow=c(2,1))
plot(data[-1,1]+0, type="l", xlim=c(1,n+t)); abline(h=nest$mu[1])
lines((vest$data[-1,]-vest$residuals)[,1]+mu1, col=2)
lines((n+1):(n+t), zt[(n+1):(n+t), 1]+mu1, lty =2)
lines((n+1):(n+t), pred$pred[,1]+mu1, lty =2, col=2)
plot(data[-1,2]+0, type="l", xlim=c(1,n+t)); abline(h=nest$mu[2])
lines((vest$data[-1,]-vest$residuals)[,2]+mu2, col=2)

mu1;mu2
