rm(list=ls())

library(openxlsx)
library(forecast)
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
library(mvtnorm)

library(ompr)
library(gurobi)
library(tidyverse)
library(ompr.roi)
library(ROI)            # load package
##remotes::install_github("FlorianSchwendinger/ROI.plugin.gurobi", INSTALL_opts="--no-multiarch", force = T)
library(roi.plugin)
ROI_installed_solvers() 
ROI_available_solvers() # check available solvers
ROI_installed_solvers() # check installed solvers


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

k <- logical()
for(j in 1:ncol(data0ar)){
  k <- c(k, adf.test(data0ar[,j])$p.value < 0.05)
}

data0arEst <- data0ar[,k]

rownames(data0arEst) <- seq(as.Date("2016/5/1"), by="month",
                    length=40)

x.orig <- data0arEst[,c(2,6,5)]
VARorder(x.orig, maxp =8)

names(x.orig) <- c("A", "B", "C")


par(mfrow=c(3,1))
for(i in 1:3){
  plot.ts(x.orig[,i],main=colnames(x.orig)[i],xlab=NULL,
          ylab='units'); grid()
}

## Unit root test
adf.test(x.orig[,1])
adf.test(x.orig[,2])
adf.test(x.orig[,3])

ccm(x.orig)
x.pacf <- pacf(x.orig, 12)
eccmz=Eccm(x.orig, maxp = 3, maxq=4, include.mean = T)
VARorder(x.orig, maxp =9)


## Vector time series model fit
m3=VAR(x.orig, p=3)
refVAR(m3, thres = 1.96)

## Diagnostic checking
m3$residuals
eccmres=Eccm(m3$residuals, maxp = 3, maxq=3)
## Vector time series model forecasting
ff <- VARpred(m3, h = 3)
ff$pred
# Upper 95%
upper<-ff$pred+1.96*ff$se
upper
# Lower 95%
lower<-ff$pred-1.96*ff$se
lower

### Simulación


rm(k)
vsim <- function(data, nobs, cnst, phi, sigma, arlags, seed = NULL){

skip <- nrow(x.orig)
  
k = nrow(sigma)
nT = nobs + skip
if(!is.null(seed)){
  set.seed(seed)
}
at = rmvnorm(nobs, rep(0, k), sigma)
nar = length(arlags)
p = 0
if (nar > 0) {
  arlags = sort(arlags)
  p = arlags[nar]
}

ist = p + 1
zt = matrix(0, nobs, k)
colnames(zt) <- names(data)
if (length(cnst) == 0) 
  cnst = rep(0, k)
zt <- as.matrix(rbind(data, zt))
for (it in (skip+1):nT) {
  tmp = matrix(at[it-skip, ], 1, k)
  if (nar > 0) {
    for (i in 1:nar) {
      idx = (i - 1) * k
      phj = phi[, (idx + 1):(idx + k)]
      ztm = matrix(zt[it - arlags[i], ], 1, k)
      tmp = tmp + ztm %*% t(phj)
    }
  }
  zt[it, ] = cnst + tmp
}
return(tail(round(zt), nobs))
}

minMax <- function(x) {
  (x - mean(x)) / (sd(x))
}

sim <- simn <- list()
for(j in 1:10000){
  cat("\r", j)
  sim[[j]] <- as.data.frame(vsim(x.orig, 6, m3$Ph0, m3$Phi, m3$Sigma, 1:3, seed = j))
simn[[j]] <- as.data.frame(matrix(c(as.matrix(as.data.frame(lapply(as.data.frame(sim[[j]]), minMax)))), nrow=1))
}

simn <- data.table::rbindlist(simn)
sim <- data.table::rbindlist(sim)

### k means

km <- kmeans(simn, centers = 10, iter.max = 10000)
library(dtwclust)

matrices <- list()
for(i in 1:1000){
  matrices[[i]] <- matrix(as.matrix(simn)[i,], nrow=6, byrow=F)
}

kmd <- tsclust(matrices, k=10)

cent <- matrix(unlist(kmd@centroids), nrow=10, byrow = T)
plot(kmd@centroids[[1]][,2], type="l", ylim=c(-2,2))
for(j in 2:10){
  lines(kmd@centroids[[j]][,2], col = j)
}

plot(km$centers[1,1:6], type="l", ylim=c(-2,2))
for(j in 2:10){
  lines(km$centers[j,1:6], col = j)
}


centroides <- list()
nn <- numeric()
for(clu in 1:10){
  idClu <- which(km$cluster == clu)
  aux <- auxVec <- list()
  for(d in 1:length(idClu)){
    aux[[d]] <- sim[[idClu[d]]]
    auxVec[[d]] <- as.data.frame(matrix(c(aux[[d]]), nrow=1))
  }
  nn[clu] <- length(aux)
  assign(paste0("DD",clu), aux)
  centroides[[clu]] <- round(as.data.frame(matrix(colMeans(data.table::rbindlist(auxVec)), nrow=1)))
  rm(aux); rm(auxVec)
}
centroides <- data.table::rbindlist(centroides)

cp1 <- as.matrix(centroides[,1:6])
cp2 <- as.matrix(centroides[,7:12])
cp3 <- as.matrix(centroides[,13:18])

dd <- array(c(cp1, cp2, cp3), dim = c(10, 6, 3))

### modelo lot sizing

oo <- c(20, 12, 1.7)
pp <- c(19.907, 12.227, 1.648)
hh <- c(2, 1.2, 0.5)
ss <- c(190, 120, 16)
CC <- c(19200, 220000, 14000)
ww <- nn/10000

model <- MIPModel() %>%
  # 1 iff i gets assigned to warehouse j
  add_variable(Q[i, t], i = 1:3, t = 1:6, type = "integer") %>%
  
  # 1 iff warehouse j is built
  add_variable(Z[i, t], i = 1:3, t = 1:6, type = "binary") %>%
  
  # 1 iff warehouse j is built
  add_variable(I[i, t, w], i = 1:3, t = 1:6, w = 1:10, type = "continuous") %>%
  
  # 1 iff warehouse j is built
  add_variable(S[i, t, w], i = 1:3, t = 1:6, w = 1:10, type = "continuous")%>%
  
  # maximize the preferences
  set_objective(sum_expr(ww[w]*oo[i] * Z[i, t] + 
                  ww[w]*pp[i] * Q[i, t] +
                  ww[w]*hh[i] * I[i, t, w] +
                  ww[w]*ss[i] * S[i, t, w], i = 1:3, t = 1:6, w = 1:10), "min") %>%
  
  # every customer needs to be assigned to a warehouse
  add_constraint(Q[1, 1] - (I[1, 1, w] - S[1, 1, w]) == dd[w,1,1], w=1:10) %>% 
  add_constraint(Q[1, 6] + (I[1, 6-1, w] - S[1, 6-1, w]) - (I[1, 6, w] - S[1, 6, w]) == dd[w, 6, 1], w=1:10) %>% 
  add_constraint(Q[1, 5] + (I[1, 5-1, w] - S[1, 5-1, w]) - (I[1, 5, w] - S[1, 5, w]) == dd[w, 5, 1], w=1:10) %>% 
  add_constraint(Q[1, 4] + (I[1, 4-1, w] - S[1, 4-1, w]) - (I[1, 4, w] - S[1, 4, w]) == dd[w, 4, 1], w=1:10) %>% 
  add_constraint(Q[1, 3] + (I[1, 3-1, w] - S[1, 3-1, w]) - (I[1, 3, w] - S[1, 3, w]) == dd[w, 3, 1], w=1:10) %>% 
  add_constraint(Q[1, 2] + (I[1, 2-1, w] - S[1, 2-1, w]) - (I[1, 2, w] - S[1, 2, w]) == dd[w, 2, 1], w=1:10) %>% 
  
  add_constraint(Q[2, 1] - (I[2, 1, w] - S[2, 1, w]) == dd[w,1,2], w=1:10) %>% 
  add_constraint(Q[2, 6] + (I[2, 6-1, w] - S[2, 6-1, w]) - (I[2, 6, w] - S[2, 6, w]) == dd[w, 6, 2], w=1:10) %>% 
  add_constraint(Q[2, 5] + (I[2, 5-1, w] - S[2, 5-1, w]) - (I[2, 5, w] - S[2, 5, w]) == dd[w, 5, 2], w=1:10) %>% 
  add_constraint(Q[2, 4] + (I[2, 4-1, w] - S[2, 4-1, w]) - (I[2, 4, w] - S[2, 4, w]) == dd[w, 4, 2], w=1:10) %>% 
  add_constraint(Q[2, 3] + (I[2, 3-1, w] - S[2, 3-1, w]) - (I[2, 3, w] - S[2, 3, w]) == dd[w, 3, 2], w=1:10) %>% 
  add_constraint(Q[2, 2] + (I[2, 2-1, w] - S[2, 2-1, w]) - (I[2, 2, w] - S[2, 2, w]) == dd[w, 2, 2], w=1:10) %>% 
  
  add_constraint(Q[3, 1] - (I[3, 1, w] - S[3, 1, w]) == dd[w,1,3], w=1:10) %>% 
  add_constraint(Q[3, 6] + (I[3, 6-1, w] - S[3, 6-1, w]) - (I[3, 6, w] - S[3, 6, w]) == dd[w, 6, 3], w=1:10) %>% 
  add_constraint(Q[3, 5] + (I[3, 5-1, w] - S[3, 5-1, w]) - (I[3, 5, w] - S[3, 5, w]) == dd[w, 5, 3], w=1:10) %>% 
  add_constraint(Q[3, 4] + (I[3, 4-1, w] - S[3, 4-1, w]) - (I[3, 4, w] - S[3, 4, w]) == dd[w, 4, 3], w=1:10) %>% 
  add_constraint(Q[3, 3] + (I[3, 3-1, w] - S[3, 3-1, w]) - (I[3, 3, w] - S[3, 3, w]) == dd[w, 3, 3], w=1:10) %>% 
  add_constraint(Q[3, 2] + (I[3, 2-1, w] - S[3, 2-1, w]) - (I[3, 2, w] - S[3, 2, w]) == dd[w, 2, 3], w=1:10) %>%
  
  add_constraint(pp[i] * Q[i, t] <= CC[i] * Z[i, t], i = 1:3, t = 1:6) %>% 
  
  add_constraint(Q[i, t] >= 0, i=1:3, t = 1:6) %>% 
  add_constraint(I[i, t, w] >= 0, i = 1:3, t = 1:6, w = 1:10) %>% 
  add_constraint(S[i, t, w] >= 0, i = 1:3, t = 1:6, w = 1:10)
  
model

result <- model %>% 
  solve_model(with_ROI(solver = "gurobi", verbose = TRUE))
result

result$solution

