rm(list=ls())

library(openxlsx)
library(forecast)

data <- read.xlsx("tasa de demanda_ejemplo (1).xlsx")

### Quitando demandas con 0
k <- logical()
for(j in 2:109){
  aux <- data[,j] == 0
  k <- c(k, sum(aux) == 0)
}

data <- data[,-1]
data0 <- data[,k]

# Manteniendo solo demandas autorregresivas
k <- logical()
for(j in 1:ncol(data0)){
  k <- c(k, sum(names(auto.arima(data0[,j])$coef) %in% "ar1") == 1)
}

data0ar <- data0[,k]

### Manteniendo variables con varianza constante

data0arc <- data0ar[,c(2,8,10,14,17)]

if("tseries" %in% rownames(installed.packages()) == FALSE)
{install.packages("tseries")}
if("MTS" %in% rownames(installed.packages()) == FALSE)
{install.packages("MTS")}
if("ppcor" %in% rownames(installed.packages()) == FALSE)
{install.packages("ppcor")}
library(tseries)
library(MTS)
library(ppcor)
library(TSA)
library(forecast)
library(ggplot2)
## Read in original data set

ms <- data0ar

#ms <- as.data.frame(read.csv("C:/Bookdata/WW2a.csv")[,-1])
rownames(ms) <- seq(as.Date("2016/5/1"), by="month",
                    length=40)
x.orig <- data0ar[,c(5,11,18)]
##Plot original vector time series
par(mfrow=c(3,1))
for(i in 1:3){
  plot.ts(x.orig[,i],main=colnames(x.orig)[i],xlab=NULL,
          ylab='units')}
## Unit root test
adf.test(x.orig[,1])
adf.test(x.orig[,2])
adf.test(x.orig[,3])

x.ccm<-ccm(x.orig, lag=15)

a <- TRUE
while(a){
  r <- sample(1:18, size = 3)
  X <- data0ar[,r]
  a <- ccm(diffM(X), level=T, output=T)$pvalue[1] < 0.05
}

## Model building
## Take difference
zt=diffM(x.orig)
## Correlation and partial correlation matrix functions
## CCM with default lag=12
zt.ccm<-ccm(zt, lag=12)


## Modeling on a revised vector series
ms <- as.data.frame(read.csv("C:/Data/Bookdata/WW2b.csv")
                    [,-1])
rownames(ms) <- seq(as.Date("2009/6/1"), by="month",
                    length=90)
x.orig <- data0arc
## Plot the new five series
par(mfrow=c(2,3))
for(i in 1:5){
  plot.ts(x.orig[,i],main=colnames(x.orig)[i],xlab=NULL,
          ylab='Mil. dollars')}
plot(seq(as.Date("2009/6/1"), by="month", length=90),x.orig
     [,1],type='l',ylab="Retail Sales",xlab="Date",ylim=c(min(x.
                                                              orig),max(x.orig)))
for(i in 2:5){
  lines(seq(as.Date("2009/6/1"), by="month", length=90),
        x.orig[,i],type='l',lty=i,lwd=2)
}
legend("topleft",legend=colnames(zt),,lty=1:5)
## Differencing the five series
ms5<-ms[,1:5]
zt=diffM(ms5)
## Plot the five differenced series
par(mfrow=c(2,3))
for(i in 1:5){
  plot.ts(zt[,i],main=colnames(x.orig)[i],xlab=NULL,
          ylab='Mil. dollars')}
# Model identication
## CCM with default

zt.pacf<-pacf(zt, 15)
eccmz=Eccm(zt)
## Vector time series model fit
m4=VAR(zt, p=3, include.mean=TRUE)
## Diagnostic checking
m4$residuals
eccmres=Eccm(m4$residuals)
## Vector time series model forecasting
ff <- VARpred(m4, h = 3)
ff$pred
# Upper 95%
upper<-ff$pred+1.96*ff$se
upper
# Lower 95%
lower<-ff$pred-1.96*ff$se
lower
q()
