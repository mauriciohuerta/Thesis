library(Rfast)
library(forecast)
library(MTS)


Sim <- Preds <- info <- list()
i <- 0
while( i < 1000){
  cat("\riteracion", i)
### Simulcaciones

 n <- 400
t <- 20

### parametros de la serie temporal

mu1 <- runif(1, 500, 1200)
mu2 <- runif(1, 500, 1200)

p11 <- runif(1, -1, 1)
p12 <- runif(1, -1, 1)
p21 <- runif(1, -1, 1)
p22 <- runif(1, -1, 1)

q11 <- runif(1, -1, 1)
q12 <- runif(1, -1, 1)
q21 <- runif(1, -1, 1)
q22 <- runif(1, -1, 1)

s11 <- runif(1, 20, 200)
s22 <- runif(1, 20, 200)
s12 <- 0

p1=matrix(c(p11,p12,p21,p22),2,2)
sig=matrix(c(s11,s12,s12,s22),2,2)
th1=matrix(c(q11,q12,q21,q22),2,2)

### Simulacion de la serie y valores reales futuros
m1=VARMAsim(n+t,arlags=c(1),malags=c(1),phi=p1,theta=th1,sigma=sig)
zt=m1$series

if(!(max(abs(zt[,1])) > 3e+5 | max(abs(zt[,2])) > 3e+5)){

#Estimacion
  
data <- zt[1:n,]

vest <- VARMACpp(data, p=1, q=1, include.mean = F)
nest <- mvnorm.mle(data)

data[,1] <- data[,1]+mu1
data[,2] <- data[,2]+mu2

# predicciones con VARMA
pred <- VARMApred(vest, h = t) 
pred$pred[,1] <- pred$pred[,1] + mu1
pred$pred[,2] <- pred$pred[,2] + mu2

i <- i+1

Sim[[i]] <- as.data.frame(t(zt))
Preds[[i]] <- as.data.frame(t(pred$pred))

info[[i]] <- list(mu=c(mu1, mu2),
                  AR = p1,
                  MA = th1,
                  Sigma = sig,
                  Vest = vest,
                  Iest = nest)
}
}

pdf("Simulaciones.pdf")
par(mfrow=c(2,1))
for(i in 1:1000){
plot(as.numeric(Sim[[i]][1,])+info[[i]]$mu[1], type="l", main=i)
plot(as.numeric(Sim[[i]][2,])+info[[i]]$mu[2], type="l")
}
dev.off()

save(Sim, file="Sim.RData")
save(Preds, file="Preds.RData")
save(info, file="info.RData")

######### Generando simulaciones futuras

rm(list=ls())

load("Sim.RData")
load("Preds.RData")
load("info.RData")

unit.sim <- function(Info.sim, Sim.sim, nT){
sigma.fit <- Info.sim$Vest$Sigma
zt <- orig <- t(Sim.sim)[1:400,]
k <- nrow(sigma.fit)
at <- rbind(matrix(0, 1, k),Info.sim$Vest$residuals)
at <- rbind(at, rmvnorm(nT, rep(0, k), sigma.fit))
zt <- as.matrix(rbind(zt, matrix(0, nT, k)))
thej <- Info.sim$Vest$Theta
phj <- Info.sim$Vest$Phi

for (it in (nrow(orig)+1):nrow(zt)){
  tmp = matrix(at[it, ], 1, k)
  atm = matrix(at[it - 1, ], 1, k)
  tmp = tmp - atm %*% t(thej)

  ztm = matrix(zt[it - 1, ], 1, k)
  tmp = tmp + ztm %*% t(phj)

  zt[it, ] = tmp
}
return(as.data.frame(matrix(as.numeric(tail(zt, nT)), nrow=1)))
}

Future.sims <- list()
for(j in 1:1000){
  sim.aux <- list()
  for(h in 1:1000){
    cat("\rSim",j,"pred",h,"        ")
  sim.aux[[h]] <- unit.sim(info[[j]], Sim[[j]], 10)
  }
  Future.sims[[j]] <- data.table::rbindlist(sim.aux)
}

save(Future.sims, file="Future.sims.RData")

##### Clustering

load("Future.sims.RData")
library(cluster)

medoids <- probs <- list()
for(k in 1:1000){
  cat("\r",k,"     ")
cl.i <- pam(Future.sims[[k]], k=20, metric = "euclidean", stand = F)
meds.i <- cl.i$medoids
prob.i <- cl.i$clusinfo[,1]/1000
medoids[[k]] <- meds.i
probs[[k]] <- prob.i
}

save(medoids, file="medoids.RData")
save(probs, file="probs.RData")

### Parametros del modelo de optimizacion

load("medoids.RData")
load("probs.RData")

info.optim <- list()
for(k in 1:1000){
  cat("\r",k,"    ")
oo <- runif(2, 250, 6000)
pp <- runif(2, 20, 30)
hh <- runif(2, 1, 10)
ss <- runif(2, 200, 300)
CC <- (info[[k]]$mu * 1.8)*pp
ww <- probs[[k]]

cp1 <- round(as.matrix(medoids[[k]][,1:10]) + info[[k]]$mu[1])
cp2 <- round(as.matrix(medoids[[k]][,11:20]) + info[[k]]$mu[2])

dd <- array(c(cp1, cp2), dim = c(20, 10, 2))

model <- MIPModel() %>%
  # 1 iff i gets assigned to warehouse j
  add_variable(Q[i, t], i = 1:2, t = 1:10, type = "integer") %>%
  
  # 1 iff warehouse j is built
  add_variable(Z[i, t], i = 1:2, t = 1:10, type = "binary") %>%
  
  # 1 iff warehouse j is built
  add_variable(I[i, t, w], i = 1:2, t = 1:10, w = 1:20, type = "continuous") %>%
  
  # 1 iff warehouse j is built
  add_variable(S[i, t, w], i = 1:2, t = 1:10, w = 1:20, type = "continuous")%>%
  
  # maximize the preferences
  set_objective(sum_expr(ww[w]*oo[i] * Z[i, t] + 
                           ww[w]*pp[i] * Q[i, t] +
                           ww[w]*hh[i] * I[i, t, w] +
                           ww[w]*ss[i] * S[i, t, w], i = 1:2, t = 1:10, w = 1:20), "min") %>%
  
  # every customer needs to be assigned to a warehouse
  add_constraint(Q[1, 1] - (I[1, 1, w] - S[1, 1, w]) == dd[w,1,1], w=1:20) %>% 
  add_constraint(Q[1, 10] + (I[1, 10-1, w] - S[1, 10-1, w]) - (I[1, 10, w] - S[1, 10, w]) == dd[w, 10, 1], w=1:20) %>%
  add_constraint(Q[1, 9] + (I[1, 9-1, w] - S[1, 9-1, w]) - (I[1, 9, w] - S[1, 9, w]) == dd[w, 9, 1], w=1:20) %>%
  add_constraint(Q[1, 8] + (I[1, 8-1, w] - S[1, 8-1, w]) - (I[1, 8, w] - S[1, 8, w]) == dd[w, 8, 1], w=1:20) %>%
  add_constraint(Q[1, 7] + (I[1, 7-1, w] - S[1, 7-1, w]) - (I[1, 7, w] - S[1, 7, w]) == dd[w, 7, 1], w=1:20) %>%
  add_constraint(Q[1, 6] + (I[1, 6-1, w] - S[1, 6-1, w]) - (I[1, 6, w] - S[1, 6, w]) == dd[w, 6, 1], w=1:20) %>% 
  add_constraint(Q[1, 5] + (I[1, 5-1, w] - S[1, 5-1, w]) - (I[1, 5, w] - S[1, 5, w]) == dd[w, 5, 1], w=1:20) %>% 
  add_constraint(Q[1, 4] + (I[1, 4-1, w] - S[1, 4-1, w]) - (I[1, 4, w] - S[1, 4, w]) == dd[w, 4, 1], w=1:20) %>% 
  add_constraint(Q[1, 3] + (I[1, 3-1, w] - S[1, 3-1, w]) - (I[1, 3, w] - S[1, 3, w]) == dd[w, 3, 1], w=1:20) %>% 
  add_constraint(Q[1, 2] + (I[1, 2-1, w] - S[1, 2-1, w]) - (I[1, 2, w] - S[1, 2, w]) == dd[w, 2, 1], w=1:20) %>% 
  
  add_constraint(Q[2, 1] - (I[2, 1, w] - S[2, 1, w]) == dd[w,1,2], w=1:20) %>% 
  add_constraint(Q[2, 10] + (I[2, 10-1, w] - S[2, 10-1, w]) - (I[2, 10, w] - S[2, 10, w]) == dd[w, 10, 2], w=1:20) %>% 
  add_constraint(Q[2, 9] + (I[2, 9-1, w] - S[2, 9-1, w]) - (I[2, 9, w] - S[2, 9, w]) == dd[w, 9, 2], w=1:20) %>% 
  add_constraint(Q[2, 8] + (I[2, 8-1, w] - S[2, 8-1, w]) - (I[2, 8, w] - S[2, 8, w]) == dd[w, 8, 2], w=1:20) %>% 
  add_constraint(Q[2, 7] + (I[2, 7-1, w] - S[2, 7-1, w]) - (I[2, 7, w] - S[2, 7, w]) == dd[w, 7, 2], w=1:20) %>% 
  add_constraint(Q[2, 6] + (I[2, 6-1, w] - S[2, 6-1, w]) - (I[2, 6, w] - S[2, 6, w]) == dd[w, 6, 2], w=1:20) %>% 
  add_constraint(Q[2, 5] + (I[2, 5-1, w] - S[2, 5-1, w]) - (I[2, 5, w] - S[2, 5, w]) == dd[w, 5, 2], w=1:20) %>% 
  add_constraint(Q[2, 4] + (I[2, 4-1, w] - S[2, 4-1, w]) - (I[2, 4, w] - S[2, 4, w]) == dd[w, 4, 2], w=1:20) %>% 
  add_constraint(Q[2, 3] + (I[2, 3-1, w] - S[2, 3-1, w]) - (I[2, 3, w] - S[2, 3, w]) == dd[w, 3, 2], w=1:20) %>% 
  add_constraint(Q[2, 2] + (I[2, 2-1, w] - S[2, 2-1, w]) - (I[2, 2, w] - S[2, 2, w]) == dd[w, 2, 2], w=1:20) %>% 
  
  add_constraint(pp[i] * Q[i, t] <= CC[i] * Z[i, t], i = 1:2, t = 1:10) %>% 
  
  add_constraint(Q[i, t] >= 0, i=1:2, t = 1:10) %>% 
  add_constraint(I[i, t, w] >= 0, i = 1:2, t = 1:10, w = 1:20) %>% 
  add_constraint(S[i, t, w] >= 0, i = 1:2, t = 1:10, w = 1:20)

model

result <- model %>% 
  solve_model(with_ROI(solver = "gurobi", verbose = TRUE))
result

result$solution

info.optim[[k]] <- list(oo = oo,
                   pp = pp,
                   hh = hh,
                   ss = ss,
                   CC = CC,
                   ww = ww,
                   dd = dd,
                   objective = result$objective_value,
                   solution = result$solution,
                   status = result$status)
}

save(info.optim, file="Info.optim.RData")

load("Info.optim.RData")

oo.s <- info.optim[[1]]$oo
pp.s = info.optim[[1]]$pp
hh.s = info.optim[[1]]$hh
ss.s = info.optim[[1]]$ss
CC.s = info.optim[[1]]$CC
ww.s = info.optim[[1]]$ww
dd.s = info.optim[[1]]$dd
objective.s = info.optim[[1]]$objective
solution.s = info.optim[[1]]$solution
status.s = info.optim[[1]]$status
ar <- as.numeric(info[[1]]$Vest$Phi)
ma <- as.numeric(info[[1]]$Vest$Theta)
si <- as.numeric(info[[1]]$Vest$Sigma)

r.ar <- as.numeric(info[[1]]$AR   )
r.ma <- as.numeric(info[[1]]$MA   )
r.si <- as.numeric(info[[1]]$Sigma)



for(k in 2:1000){
  oo.s <- c(oo.s, info.optim[[k]]$oo)
  pp.s <- c(pp.s, info.optim[[k]]$pp)
  hh.s <- c(hh.s, info.optim[[k]]$hh)
  ss.s <- c(ss.s, info.optim[[k]]$ss)
  CC.s <- c(CC.s, info.optim[[k]]$CC)
  ww.s <- c(ww.s, info.optim[[k]]$ww)
  objective.s <- c(objective.s, info.optim[[k]]$objective)
  solution.s <- c(solution.s, info.optim[[k]]$solution)
  ar <- c(ar, as.numeric(info[[k]]$Vest$Phi))
  ma <- c(ma, as.numeric(info[[k]]$Vest$Theta))
  si <- c(si, as.numeric(info[[k]]$Vest$Sigma))
  
  r.ar <- c(r.ar, as.numeric(info[[k]]$AR   ))
  r.ma <- c(r.ma, as.numeric(info[[k]]$MA   ))
  r.si <- c(r.si, as.numeric(info[[k]]$Sigma))
  }

oo.s <- matrix(oo.s, ncol=2, byrow = T)
pp.s <- matrix(pp.s, ncol=2, byrow = T)
hh.s <- matrix(hh.s, ncol=2, byrow = T)
ss.s <- matrix(ss.s, ncol=2, byrow = T)
CC.s <- matrix(CC.s, ncol=2, byrow = T)
ww.s <- matrix(ww.s, ncol=20, byrow = T)
solution.s <- matrix(solution.s, ncol=840, byrow = T)
ar <- matrix(ar, ncol=4, byrow = T)
ma <- matrix(ma, ncol=4, byrow = T)
si <- matrix(si, ncol=4, byrow = T)
r.ar <- matrix(r.ar, ncol=4, byrow = T)
r.ma <- matrix(r.ma, ncol=4, byrow = T)
r.si <- matrix(r.si, ncol=4, byrow = T)


summary(lm(objective.s~oo.s[,1]+oo.s[,2]))
summary(lm(objective.s~pp.s[,1]+pp.s[,2]))
summary(lm(objective.s~hh.s[,1]+hh.s[,2]))
summary(lm(objective.s~ss.s[,1]+ss.s[,2]))
summary(lm(objective.s~CC.s[,1]+CC.s[,2]))
summary(lm(objective.s~ww.s))
summary(lm(objective.s~abs(ar)))
summary(lm(objective.s~abs(ma)))
summary(lm(objective.s~si[,-3]))

summary(lm(objective.s~abs(r.ar)))
summary(lm(objective.s~abs(r.ma)))
summary(lm(objective.s~r.si[,c(1,4)]))

summary(lm(objective.s~CC.s+oo.s+pp.s+hh.s+ss.s+abs(r.ar)+abs(r.ma)))

par(mfrow=c(1,2))
plot(oo.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(oo.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(1,2))
plot(pp.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(pp.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(1,2))
plot(hh.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(hh.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(1,2))
plot(ss.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(ss.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(1,2))
plot(CC.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(CC.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(2,2))
plot((ar[,1]), objective.s/1000, pch=16, cex=0.7); grid()
plot((ar[,3]), objective.s/1000, pch=16, cex=0.7); grid()
plot((ar[,2]), objective.s/1000, pch=16, cex=0.7); grid()
plot((ar[,4]), objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(2,2))
plot((r.ar[,1]), objective.s/1000, pch=16, cex=0.7); grid()
plot((r.ar[,3]), objective.s/1000, pch=16, cex=0.7); grid()
plot((r.ar[,2]), objective.s/1000, pch=16, cex=0.7); grid()
plot((r.ar[,4]), objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(2,2))
plot(abs(ma[,1]), objective.s/1000, pch=16, cex=0.7); grid()
plot(abs(ma[,3]), objective.s/1000, pch=16, cex=0.7); grid()
plot(abs(ma[,2]), objective.s/1000, pch=16, cex=0.7); grid()
plot(abs(ma[,4]), objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(1,2))
plot(CC.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(CC.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(1,2))
plot(CC.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(CC.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

###### Costos marginales

marCost <- function(infor){
sol <- infor$solution
costos <- c(rep(infor$pp, each=10),
            rep(infor$oo, each = 10),
            rep(infor$hh, each=10*20),
            rep(infor$ss, each=10*20))
pesos <- c(rep(1, 10*2*2), 
           rep(infor$ww, times = 10*2*2))

Ccost <- pesos*costos*sol
complete.marginal.cost <- matrix(Ccost, nrow = 1)
colnames(complete.marginal.cost) <- names(sol)

E.cost <- Ccost[1:40]
for(i in seq(40,820, by=20)){
  E.cost <- c(E.cost, sum(Ccost[(i+1):(i+20)]))
}
                
namE <- c(rep(paste0("Q.", 1:2), each=10),
          rep(paste0("Z.", 1:2), each=10),
          rep(paste0("I.", 1:2), each=10),
          rep(paste0("S.", 1:2), each=10))
namE <- paste0(namE, "_",formatC(1:10, width=2, flag = "0"))

E.cost <- matrix(E.cost, nrow = 1)
colnames(E.cost) <- namE
return(list(Complete = complete.marginal.cost,
            Expected = E.cost))
}

EC <- marCost(info.optim[[1]])$Expected
for(k in 2:1000){
  cat("\r",k,"    ")
  EC <- rbind(EC, marCost(info.optim[[k]])$Expected)
}

EC1 <- EC[,seq(1,80, by=10)]
EC1 <- rowSums(EC1)

summary(lm(EC1~CC.s+oo.s+pp.s+hh.s+ss.s+abs(ar)+abs(ma)))


summary(lm(EC[,71]~abs(ar)))
plot(CC.s[,2], EC1)


##### Replicacion sin temporalidad



######### Generando simulaciones futuras

rm(list=ls())

load("Sim.RData")
load("Preds.RData")
load("info.RData")

unit.sim <- function(Info.sim, Sim.sim, nT){
  sigma.fit <- Info.sim$Iest$sigma
  mu.fit <- Info.sim$Iest$mu + info[[1]]$mu
  zt <- orig <- t(Sim.sim)[1:400,]
  k <- nrow(sigma.fit)
  zt <- as.matrix(rbind(zt, Rfast::rmvnorm(nT, mu.fit, sigma = sigma.fit)))
  return(as.data.frame(matrix(as.numeric(tail(zt, nT)), nrow=1)))
}

Future.sims <- list()
for(j in 1:1000){
  sim.aux <- list()
  for(h in 1:1000){
    cat("\rSim",j,"pred",h,"        ")
    sim.aux[[h]] <- unit.sim(info[[j]], Sim[[j]], 10)
  }
  Future.sims[[j]] <- data.table::rbindlist(sim.aux)
}

save(Future.sims, file="Ind/Future.sims.RData")
load("Ind/Future.sims.RData")


##### Clustering

library(cluster)

medoids <- probs <- list()
for(k in 1:1000){
  cat("\r",k,"     ")
  cl.i <- pam(Future.sims[[k]], k=20, metric = "euclidean", stand = F)
  meds.i <- cl.i$medoids
  prob.i <- cl.i$clusinfo[,1]/1000
  medoids[[k]] <- meds.i
  probs[[k]] <- prob.i
}

save(medoids, file="Ind/medoids.RData")
save(probs, file="Ind/probs.RData")

### Parametros del modelo de optimizacion

info.optim <- list()
for(k in 1:1000){
  cat("\r",k,"    ")
  oo <- runif(2, 250, 6000)
  pp <- runif(2, 20, 30)
  hh <- runif(2, 1, 10)
  ss <- runif(2, 200, 300)
  CC <- (info[[k]]$mu * 1.8)*pp
  ww <- probs[[k]]
  
  cp1 <- round(as.matrix(medoids[[k]][,1:10]) + info[[k]]$mu[1])
  cp2 <- round(as.matrix(medoids[[k]][,11:20]) + info[[k]]$mu[2])
  
  dd <- array(c(cp1, cp2), dim = c(20, 10, 2))
  
  model <- MIPModel() %>%
    # 1 iff i gets assigned to warehouse j
    add_variable(Q[i, t], i = 1:2, t = 1:10, type = "integer") %>%
    
    # 1 iff warehouse j is built
    add_variable(Z[i, t], i = 1:2, t = 1:10, type = "binary") %>%
    
    # 1 iff warehouse j is built
    add_variable(I[i, t, w], i = 1:2, t = 1:10, w = 1:20, type = "continuous") %>%
    
    # 1 iff warehouse j is built
    add_variable(S[i, t, w], i = 1:2, t = 1:10, w = 1:20, type = "continuous")%>%
    
    # maximize the preferences
    set_objective(sum_expr(ww[w]*oo[i] * Z[i, t] + 
                             ww[w]*pp[i] * Q[i, t] +
                             ww[w]*hh[i] * I[i, t, w] +
                             ww[w]*ss[i] * S[i, t, w], i = 1:2, t = 1:10, w = 1:20), "min") %>%
    
    # every customer needs to be assigned to a warehouse
    add_constraint(Q[1, 1] - (I[1, 1, w] - S[1, 1, w]) == dd[w,1,1], w=1:20) %>% 
    add_constraint(Q[1, 10] + (I[1, 10-1, w] - S[1, 10-1, w]) - (I[1, 10, w] - S[1, 10, w]) == dd[w, 10, 1], w=1:20) %>%
    add_constraint(Q[1, 9] + (I[1, 9-1, w] - S[1, 9-1, w]) - (I[1, 9, w] - S[1, 9, w]) == dd[w, 9, 1], w=1:20) %>%
    add_constraint(Q[1, 8] + (I[1, 8-1, w] - S[1, 8-1, w]) - (I[1, 8, w] - S[1, 8, w]) == dd[w, 8, 1], w=1:20) %>%
    add_constraint(Q[1, 7] + (I[1, 7-1, w] - S[1, 7-1, w]) - (I[1, 7, w] - S[1, 7, w]) == dd[w, 7, 1], w=1:20) %>%
    add_constraint(Q[1, 6] + (I[1, 6-1, w] - S[1, 6-1, w]) - (I[1, 6, w] - S[1, 6, w]) == dd[w, 6, 1], w=1:20) %>% 
    add_constraint(Q[1, 5] + (I[1, 5-1, w] - S[1, 5-1, w]) - (I[1, 5, w] - S[1, 5, w]) == dd[w, 5, 1], w=1:20) %>% 
    add_constraint(Q[1, 4] + (I[1, 4-1, w] - S[1, 4-1, w]) - (I[1, 4, w] - S[1, 4, w]) == dd[w, 4, 1], w=1:20) %>% 
    add_constraint(Q[1, 3] + (I[1, 3-1, w] - S[1, 3-1, w]) - (I[1, 3, w] - S[1, 3, w]) == dd[w, 3, 1], w=1:20) %>% 
    add_constraint(Q[1, 2] + (I[1, 2-1, w] - S[1, 2-1, w]) - (I[1, 2, w] - S[1, 2, w]) == dd[w, 2, 1], w=1:20) %>% 
    
    add_constraint(Q[2, 1] - (I[2, 1, w] - S[2, 1, w]) == dd[w,1,2], w=1:20) %>% 
    add_constraint(Q[2, 10] + (I[2, 10-1, w] - S[2, 10-1, w]) - (I[2, 10, w] - S[2, 10, w]) == dd[w, 10, 2], w=1:20) %>% 
    add_constraint(Q[2, 9] + (I[2, 9-1, w] - S[2, 9-1, w]) - (I[2, 9, w] - S[2, 9, w]) == dd[w, 9, 2], w=1:20) %>% 
    add_constraint(Q[2, 8] + (I[2, 8-1, w] - S[2, 8-1, w]) - (I[2, 8, w] - S[2, 8, w]) == dd[w, 8, 2], w=1:20) %>% 
    add_constraint(Q[2, 7] + (I[2, 7-1, w] - S[2, 7-1, w]) - (I[2, 7, w] - S[2, 7, w]) == dd[w, 7, 2], w=1:20) %>% 
    add_constraint(Q[2, 6] + (I[2, 6-1, w] - S[2, 6-1, w]) - (I[2, 6, w] - S[2, 6, w]) == dd[w, 6, 2], w=1:20) %>% 
    add_constraint(Q[2, 5] + (I[2, 5-1, w] - S[2, 5-1, w]) - (I[2, 5, w] - S[2, 5, w]) == dd[w, 5, 2], w=1:20) %>% 
    add_constraint(Q[2, 4] + (I[2, 4-1, w] - S[2, 4-1, w]) - (I[2, 4, w] - S[2, 4, w]) == dd[w, 4, 2], w=1:20) %>% 
    add_constraint(Q[2, 3] + (I[2, 3-1, w] - S[2, 3-1, w]) - (I[2, 3, w] - S[2, 3, w]) == dd[w, 3, 2], w=1:20) %>% 
    add_constraint(Q[2, 2] + (I[2, 2-1, w] - S[2, 2-1, w]) - (I[2, 2, w] - S[2, 2, w]) == dd[w, 2, 2], w=1:20) %>% 
    
    add_constraint(pp[i] * Q[i, t] <= CC[i] * Z[i, t], i = 1:2, t = 1:10) %>% 
    
    add_constraint(Q[i, t] >= 0, i=1:2, t = 1:10) %>% 
    add_constraint(I[i, t, w] >= 0, i = 1:2, t = 1:10, w = 1:20) %>% 
    add_constraint(S[i, t, w] >= 0, i = 1:2, t = 1:10, w = 1:20)
  
  model
  
  result <- model %>% 
    solve_model(with_ROI(solver = "gurobi", verbose = TRUE))
  result
  
  result$solution
  
  info.optim[[k]] <- list(oo = oo,
                          pp = pp,
                          hh = hh,
                          ss = ss,
                          CC = CC,
                          ww = ww,
                          dd = dd,
                          objective = result$objective_value,
                          solution = result$solution,
                          status = result$status)
}

save(info.optim, file="Ind/Info.optim.RData")

oo.s <- info.optim[[1]]$oo
pp.s = info.optim[[1]]$pp
hh.s = info.optim[[1]]$hh
ss.s = info.optim[[1]]$ss
CC.s = info.optim[[1]]$CC
ww.s = info.optim[[1]]$ww
dd.s = info.optim[[1]]$dd
objective.s = info.optim[[1]]$objective
solution.s = info.optim[[1]]$solution
status.s = info.optim[[1]]$status
ar <- as.numeric(info[[1]]$Vest$Phi)
ma <- as.numeric(info[[1]]$Vest$Theta)
si <- as.numeric(info[[1]]$Vest$Sigma)

r.ar <- as.numeric(info[[1]]$AR   )
r.ma <- as.numeric(info[[1]]$MA   )
r.si <- as.numeric(info[[1]]$Sigma)



for(k in 2:1000){
  oo.s <- c(oo.s, info.optim[[k]]$oo)
  pp.s <- c(pp.s, info.optim[[k]]$pp)
  hh.s <- c(hh.s, info.optim[[k]]$hh)
  ss.s <- c(ss.s, info.optim[[k]]$ss)
  CC.s <- c(CC.s, info.optim[[k]]$CC)
  ww.s <- c(ww.s, info.optim[[k]]$ww)
  objective.s <- c(objective.s, info.optim[[k]]$objective)
  solution.s <- c(solution.s, info.optim[[k]]$solution)
  ar <- c(ar, as.numeric(info[[k]]$Vest$Phi))
  ma <- c(ma, as.numeric(info[[k]]$Vest$Theta))
  si <- c(si, as.numeric(info[[k]]$Vest$Sigma))
  
  r.ar <- c(r.ar, as.numeric(info[[k]]$AR   ))
  r.ma <- c(r.ma, as.numeric(info[[k]]$MA   ))
  r.si <- c(r.si, as.numeric(info[[k]]$Sigma))
}

oo.s <- matrix(oo.s, ncol=2, byrow = T)
pp.s <- matrix(pp.s, ncol=2, byrow = T)
hh.s <- matrix(hh.s, ncol=2, byrow = T)
ss.s <- matrix(ss.s, ncol=2, byrow = T)
CC.s <- matrix(CC.s, ncol=2, byrow = T)
ww.s <- matrix(ww.s, ncol=20, byrow = T)
solution.s <- matrix(solution.s, ncol=840, byrow = T)
ar <- matrix(ar, ncol=4, byrow = T)
ma <- matrix(ma, ncol=4, byrow = T)
si <- matrix(si, ncol=4, byrow = T)
r.ar <- matrix(r.ar, ncol=4, byrow = T)
r.ma <- matrix(r.ma, ncol=4, byrow = T)
r.si <- matrix(r.si, ncol=4, byrow = T)


summary(lm(objective.s~oo.s[,1]+oo.s[,2]))
summary(lm(objective.s~pp.s[,1]+pp.s[,2]))
summary(lm(objective.s~hh.s[,1]+hh.s[,2]))
summary(lm(objective.s~ss.s[,1]+ss.s[,2]))
summary(lm(objective.s~CC.s[,1]+CC.s[,2]))
summary(lm(objective.s~ww.s))
summary(lm(objective.s~abs(ar)))
summary(lm(objective.s~abs(ma)))
summary(lm(objective.s~si))

summary(lm(objective.s~abs(r.ar)))
summary(lm(objective.s~abs(r.ma)))
summary(lm(objective.s~r.si[,c(1,4)]))

summary(lm(objective.s~CC.s+oo.s+pp.s+hh.s+ss.s+abs(ar)+abs(ma)))

par(mfrow=c(1,2))
plot(oo.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(oo.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(1,2))
plot(pp.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(pp.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(1,2))
plot(hh.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(hh.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(1,2))
plot(ss.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(ss.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(1,2))
plot(CC.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(CC.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(2,2))
plot(abs(ar[,1]), objective.s/1000, pch=16, cex=0.7); grid()
plot(abs(ar[,3]), objective.s/1000, pch=16, cex=0.7); grid()
plot(abs(ar[,2]), objective.s/1000, pch=16, cex=0.7); grid()
plot(abs(ar[,4]), objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(2,2))
plot(abs(r.ar[,1]), objective.s/1000, pch=16, cex=0.7); grid()
plot(abs(r.ar[,3]), objective.s/1000, pch=16, cex=0.7); grid()
plot(abs(r.ar[,2]), objective.s/1000, pch=16, cex=0.7); grid()
plot(abs(r.ar[,4]), objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(2,2))
plot(abs(ma[,1]), objective.s/1000, pch=16, cex=0.7); grid()
plot(abs(ma[,3]), objective.s/1000, pch=16, cex=0.7); grid()
plot(abs(ma[,2]), objective.s/1000, pch=16, cex=0.7); grid()
plot(abs(ma[,4]), objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(1,2))
plot(CC.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(CC.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

par(mfrow=c(1,2))
plot(CC.s[,1], objective.s/1000, pch=16, cex=0.7); grid()
plot(CC.s[,2], objective.s/1000, pch=16, cex=0.7); grid()

###### Costos marginales

marCost <- function(infor){
  sol <- infor$solution
  costos <- c(rep(infor$pp, each=10),
              rep(infor$oo, each = 10),
              rep(infor$hh, each=10*20),
              rep(infor$ss, each=10*20))
  pesos <- c(rep(1, 10*2*2), 
             rep(infor$ww, times = 10*2*2))
  
  Ccost <- pesos*costos*sol
  complete.marginal.cost <- matrix(Ccost, nrow = 1)
  colnames(complete.marginal.cost) <- names(sol)
  
  E.cost <- Ccost[1:40]
  for(i in seq(40,820, by=20)){
    E.cost <- c(E.cost, sum(Ccost[(i+1):(i+20)]))
  }
  
  namE <- c(rep(paste0("Q.", 1:2), each=10),
            rep(paste0("Z.", 1:2), each=10),
            rep(paste0("I.", 1:2), each=10),
            rep(paste0("S.", 1:2), each=10))
  namE <- paste0(namE, "_",formatC(1:10, width=2, flag = "0"))
  
  E.cost <- matrix(E.cost, nrow = 1)
  colnames(E.cost) <- namE
  return(list(Complete = complete.marginal.cost,
              Expected = E.cost))
}

EC <- marCost(info.optim[[1]])$Expected
for(k in 2:1000){
  cat("\r",k,"    ")
  EC <- rbind(EC, marCost(info.optim[[k]])$Expected)
}

EC1 <- EC[,seq(1,80, by=10)]
EC1 <- rowSums(EC1)

summary(lm(EC1~CC.s+oo.s+pp.s+hh.s+ss.s+abs(r.ar)+abs(r.ma)))


summary(lm(EC[,2]~abs(ar)))
plot(CC.s[,2], EC1)
