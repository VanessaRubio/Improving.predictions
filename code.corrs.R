rm(list = ls())

packages = c("data.table","lmerTest","lme4","MuMIn","bbmle","psych")
lapply(packages,require,character.only=T)

initial = as.data.frame(survival)

#Log transformation
initial$P = log(initial$P)  #log of phosphorus concentration
initial$N = log(initial$N)  #log of nitrogen concentration
initial$la = log(initial$la)  #log of leaf_area concentration
initial$h = log(initial$h)  #log of max_height concentration
initial$seed = log(initial$seed)  #log of seedmass concentration
initial$sla = log(initial$sla)  #log of specific_leaf_area concentration
initial$vla = log(initial$vla) #log of vein length per area concentration

species=match(unique(initial$species),initial$species) 
initial.cor=initial[species,]
initial.cor = initial.cor[-c(27),]
rownames(initial.cor) = 1:dim(initial.cor)[1]

initial.cor$dbh = NULL
initial.cor$TNCI = NULL
initial.cor$RGR = NULL
data = initial.cor[7:16]

c13.cor <- sapply(colnames(data), function(trait.name){
	c13 <- data[["c13"]]
	trait <- data[[trait.name]]
	
	correlation.test <- cor.test(c13, trait)
	correlation <- as.vector(correlation.test$estimate)
	p.value <- as.vector(correlation.test$p.value)
	n <- as.numeric(correlation.test$parameter)
	out.mat <- c(correlation, p.value, n)
	names(out.mat) <- c("correlation", "p.value", "n")
	return(out.mat)
})
vla.cor <- sapply(colnames(data), function(trait.name){
  vla <- data[["vla"]]
  trait <- data[[trait.name]]
  
  correlation.test <- cor.test(vla, trait)
  correlation <- as.vector(correlation.test$estimate)
  p.value <- as.vector(correlation.test$p.value)
  n <- as.numeric(correlation.test$parameter)
  out.mat <- c(correlation, p.value, n)
  names(out.mat) <- c("correlation", "p.value","n")
  return(out.mat)
})
P.cor <- sapply(colnames(data), function(trait.name){
  P <- data[["P"]]
  trait <- data[[trait.name]]
  
  correlation.test <- cor.test(P, trait)
  correlation <- as.vector(correlation.test$estimate)
  p.value <- as.vector(correlation.test$p.value)
  n <- as.numeric(correlation.test$parameter)
  out.mat <- c(correlation, p.value, n)
  names(out.mat) <- c("correlation", "p.value","n")
  return(out.mat)
})
N.cor <- sapply(colnames(data), function(trait.name){
  N <- data[["N"]]
  trait <- data[[trait.name]]
  
  correlation.test <- cor.test(N, trait)
  correlation <- as.vector(correlation.test$estimate)
  p.value <- as.vector(correlation.test$p.value)
  n <- as.numeric(correlation.test$parameter)
  out.mat <- c(correlation, p.value, n)
  names(out.mat) <- c("correlation", "p.value","n")
  return(out.mat)
})
C.cor <- sapply(colnames(data), function(trait.name){
  C <- data[["C"]]
  trait <- data[[trait.name]]
  
  correlation.test <- cor.test(C, trait)
  correlation <- as.vector(correlation.test$estimate)
  p.value <- as.vector(correlation.test$p.value)
  n <- as.numeric(correlation.test$parameter)
  out.mat <- c(correlation, p.value, n)
  names(out.mat) <- c("correlation", "p.value","n")
  return(out.mat)
})
wsg.cor <- sapply(colnames(data), function(trait.name){
  wsg <- data[["wsg"]]
  trait <- data[[trait.name]]
  
  correlation.test <- cor.test(wsg, trait)
  correlation <- as.vector(correlation.test$estimate)
  p.value <- as.vector(correlation.test$p.value)
  n <- as.numeric(correlation.test$parameter)
  out.mat <- c(correlation, p.value, n)
  names(out.mat) <- c("correlation", "p.value","n")
  return(out.mat)
})
la.cor <- sapply(colnames(data), function(trait.name){
  la <- data[["la"]]
  trait <- data[[trait.name]]
  
  correlation.test <- cor.test(la, trait)
  correlation <- as.vector(correlation.test$estimate)
  p.value <- as.vector(correlation.test$p.value)
  n <- as.numeric(correlation.test$parameter)
  out.mat <- c(correlation, p.value, n)
  names(out.mat) <- c("correlation", "p.value","n")
  return(out.mat)
})
h.cor <- sapply(colnames(data), function(trait.name){
  h <- data[["h"]]
  trait <- data[[trait.name]]
  
  correlation.test <- cor.test(h, trait)
  correlation <- as.vector(correlation.test$estimate)
  p.value <- as.vector(correlation.test$p.value)
  n <- as.numeric(correlation.test$parameter)
  out.mat <- c(correlation, p.value, n)
  names(out.mat) <- c("correlation", "p.value","n")
  return(out.mat)
})
seed.cor <- sapply(colnames(data), function(trait.name){
  seed <- data[["seed"]]
  trait <- data[[trait.name]]
  
  correlation.test <- cor.test(seed, trait)
  correlation <- as.vector(correlation.test$estimate)
  p.value <- as.vector(correlation.test$p.value)
  n <- as.numeric(correlation.test$parameter)
  out.mat <- c(correlation, p.value, n)
  names(out.mat) <- c("correlation", "p.value","n")
  return(out.mat)
})
sla.cor <- sapply(colnames(data), function(trait.name){
  sla <- data[["sla"]]
  trait <- data[[trait.name]]
  
  correlation.test <- cor.test(sla, trait)
  correlation <- as.vector(correlation.test$estimate)
  p.value <- as.vector(correlation.test$p.value)
  n <- as.numeric(correlation.test$parameter)
  out.mat <- c(correlation, p.value, n)
  names(out.mat) <- c("correlation", "p.value","n")
  return(out.mat)
})

# significant correlations between # c13 and P , and , vla and C
##
par(mfrow = c(1,3))
plot(cex.lab = 1.2, cex.axis=1.2,col = adjustcolor(9,alpha.f=0.45),cex=2.5,initial.cor$c13, initial.cor$P, ylab = "log (leaf phosphorus concentration (% P of dry mass))", xlab = "Leaf carbon isotope composition (%o)",pch = 20)
legend(cex =0.9,-27.2,-3, inset = c(0, 1), c("r = 0.27","p-value = 0.004", "n = 105"),bty = "n", xpd = T, plot = T)
m1 = lm(initial.cor$P~initial.cor$c13)
abline(m1)

plot(cex.lab = 1.2, cex.axis=1.2,col = adjustcolor(9,alpha.f=0.45),cex=2.5,initial.cor$c13, initial.cor$wsg, ylab = "Wood density", xlab = "Leaf carbon isotope composition (%o)",pch = 20)
legend(-27.8,0.99, inset = c(0, 1), c("r = -0.20","p-value = 0.03", "n = 105"),bty = "n", xpd = T, plot = T)
m2 = lm(initial.cor$wsg~initial.cor$c13)
abline(m2)

plot(cex.lab = 1.2, cex.axis=1.2,col = adjustcolor(9,alpha.f=0.45),cex=2.5,initial.cor$vla, initial.cor$C, ylab = "Leaf carbon concentration (% C of dry mass)", xlab = "log(vein length per unit area (mm mm-2))",pch = 20)
legend(4.2,41.3, inset = c(0, 1), c("r = 0.42","p-value <0.001", "n = 60"),bty = "n", xpd = T, plot = T)
m3 = lm(initial.cor$C~initial.cor$vla)
abline(m3)
