rm(list = ls())

packages = c("data.table","lmerTest","lme4","MuMIn","bbmle")
lapply(packages,require,character.only=T)

initial = as.data.frame(survival)
survival_transf = subset(initial, !is.na(vla))

#Log transformation
survival_transf$P = log(survival_transf$P)  #log of phosphorus concentration
survival_transf$N = log(survival_transf$N)  #log of nitrogen concentration
survival_transf$la = log(survival_transf$la)  #log of leaf_area concentration
survival_transf$h = log(survival_transf$h)  #log of max_height concentration
survival_transf$seed = log(survival_transf$seed)  #log of seedmass concentration
survival_transf$sla = log(survival_transf$sla)  #log of specific_leaf_area concentration

rownames(survival_transf) = 1:dim(survival_transf)[1]

survival_transf$tag = NULL
survival_transf$stemtag = NULL
survival_transf$lxm = NULL
survival_transf$lym = NULL
fo = match(unique(survival_transf$species),survival_transf$species)
forpca = survival_transf[fo,]
all.pca = forpca[8:17] #all 60 species
pca = prcomp(all.pca,retx=TRUE,scale.=T)

sd = pca$sdev
loadings = data.frame(pca$rotation)
scores = pca$x 
var = sd^2
var.percent = var/sum(var) * 100

#Broken stick criterion
col = 10   
barplot(var.percent, xlab="PC", ylab="Percent Variance", names.arg=1:length(var.percent), las=1, ylim=c(0,max(var.percent)), col="gray")
abline(h=1/col*100, col="red")

rownames(scores)=forpca[,1]

### plots
par(mfrow = c(1,3), mai = c(0.4,0.65,0.4,0.6))
biplot(pca_incl,choices = c(1,2) ,scale = T, cex = 2, cex.axis = 2.5, cex.lab = 2.5, xlabs=rep("", dim(pca_incl$x)[1]), col = c("black"), xlim = c(-0.32,0.32), ylim = c(-0.32,0.32),xlab = "PC1 (27.8%)", ylab = "PC2 (18.4%)")
biplot(pca_incl,choices = c(2,3) ,scale = T, cex = 2, cex.axis = 2.5, cex.lab = 2.5, xlabs=rep("", dim(pca_incl$x)[1]), col = c("black"), xlim = c(-0.32,0.32), ylim = c(-0.32,0.32),xlab = "PC2 (18.4%)", ylab = "PC3 (13.9%)")
biplot(pca_incl,choices = c(1,3) ,scale = T, cex = 2, cex.axis = 2.5, cex.lab = 2.5, xlabs=rep("", dim(pca_incl$x)[1]), col = c("black"), xlim = c(-0.32,0.32), ylim = c(-0.32,0.32),xlab = "PC1 (27.8%)", ylab = "PC3 (13.9%)")


