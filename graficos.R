### Resultados de la simulacion

rm(list=ls())

load("Sim.RData")
load("Preds.RData")
load("info.RData")

load("medoids.RData")
load("probs.RData")

load("Info.optim.RData")

#### Analisis

oo.s <- info.optim[[1]]$oo
pp.s = info.optim[[1]]$pp
hh.s = info.optim[[1]]$hh
ss.s = info.optim[[1]]$ss
CC.s = info.optim[[1]]$CC
ww.s = info.optim[[1]]$ww
dd.s = info[[1]]$mu
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
  dd.s <- c(dd.s, info[[k]]$mu)
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
dd.s <- matrix(dd.s, ncol=2, byrow = T)
CC.s <- matrix(CC.s, ncol=2, byrow = T)
ww.s <- matrix(ww.s, ncol=20, byrow = T)
solution.s <- matrix(solution.s, ncol=840, byrow = T)
ar <- matrix(ar, ncol=4, byrow = T)
ma <- matrix(ma, ncol=4, byrow = T)
si <- matrix(si, ncol=4, byrow = T)
r.ar <- matrix(r.ar, ncol=4, byrow = T)
r.ma <- matrix(r.ma, ncol=4, byrow = T)
r.si <- matrix(r.si, ncol=4, byrow = T)

ar <- -ar
summary(lm(objective.s~0+oo.s+pp.s+hh.s+ss.s+dd.s+CC.s+ar+ma+si[,-3]))

datos_resultados <- cbind(objective.s,oo.s,pp.s,hh.s,ss.s,dd.s,CC.s,ar,ma,si[,-2])
colnames(datos_resultados) <- c(
  "Expected cost",
  "ordering cost 1",
  "ordering cost 2",
  "unit purchase cost 1",
  "unit purchase cost 2",
  "holding cost 1",
  "holding cost 2",
  "shortage cost 1",
  "shortage cost 2",
  "mu 1",
  "mu 2",
  "capacity 1",
  "capacity 2",
  "Phi(1,1)",
  "Phi(1,2)",
  "Phi(2,1)",
  "Phi(2,2)",
  "Psi(1,1)",
  "Psi(1,2)",
  "Psi(2,1)",
  "Psi(2,2)",
  "Sigma(1,1)",
  "Sigma(1,2)",
  "Sigma(2,2)"
)
corrplot::corrplot(cor(datos_resultados))
corrplot::corrplot(cor(datos_resultados), p.mat=psych::corr.test(datos_resultados, use="pairwise", insig = "p-value")$p.value)

library(PerformanceAnalytics)
chart.Correlation(datos)


for(i in 2:ncol(datos_resultados)){
  plot(datos_resultados[,i], datos_resultados[,1], pch=16, col="steelblue", cex=0.6,
       xlab=colnames(datos_resultados)[i], ylab="expected costs")
}

### Analisis con filtro

library(tidyverse)
datos <- as.data.frame(cbind(datos_resultados, r.ar, r.ma, r.si[,-c(2,3)]))
names(datos) <- c(colnames(datos_resultados), 
                           "Ar11", "AR12", "AR21", "AR22",
                           "MA11", "MA12", "MA21", "MA22",
                           "SI11", "SI22")

datos <- datos %>% 
  mutate(r = `Sigma(1,2)`/(sqrt(`Sigma(1,1)`*`Sigma(2,2)`)),
         SI12 = r*sqrt(SI11*SI22)) %>% 
  filter(`Expected cost` < 800000,
         `Sigma(1,1)` < 1000,
         `Sigma(1,2)` < 1000,
         `Sigma(2,2)` < 1000)

plot(1/(datos$r), datos$`Expected cost`)

hist(datos$r, freq = F)
curve(dPE2(x, mu=-0.01006, sigma=.09, nu=.9), add=T)

xa <- pPE2(datos$r, mu=-0.01006, sigma=.09, nu=.9)
xb <- qunif(xa, -1, 1)

plot(abs(xb)^2, datos$`Expected cost`)

summary(lm(`Expected cost`~., data = datos[,1:24]))

####

pdf("HistTC.pdf", height = 5, width = 6)
options(scipen = 999999)
hist(datos$`Expected cost`/1000, freq=F, col="white", ylim=c(0,.005), main="",
     xlab="expected total cost (in thousands)", ylab="density")

mu.est <- mean(datos$`Expected cost`)
sd.est <- sd(datos$`Expected cost`)

curve(dnorm(x, mu.est/1000, sd.est/1000), add=T, col=2)

ks.test(datos$`Expected cost`, "pnorm", mean=mu.est, sd=sd.est)
grid()
text(200, 0.004, labels="KS test for normality\np-value: 0.7425", pos=4)
dev.off()
