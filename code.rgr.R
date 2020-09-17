### GROWTH ####

rm(list = ls())
set.seed(NULL)
packages = c("lmerTest","lme4","MuMIn","qrmix","data.table","parallel","aod", "DMwR")
lapply(packages,require,character.only=T)

growth = fread("~/Data.txt")

initial = as.data.frame(growth)
growth_transf= subset(initial, !is.na(ca2)) #deletes the species that do not have canopy values
growth_transf = subset(growth_transf, !is.na(RGR)) #deletes the NAs RGR of new individuals in census 2011
growth_transf = subset(growth_transf, RGR>=0) 
growth_transf$RGR = growth_transf$RGR+1
growth_transf = subset(growth_transf,TNCI >0 )

growth_transf$log_ca2 = log(growth_transf$ca2)
growth_transf$log_ca2 = scale(growth_transf$log_ca2)   #crown area

growth_transf$mp = log(growth_transf$ca2*growth_transf$lma)
growth_transf$mp = scale(growth_transf$mp)

#Log transformation and/or transform
growth_transf$dbh = scale(log(growth_transf$dbh))
growth_transf$dbh2 = scale(log(growth_transf$dbh2))
growth_transf$P = scale(log(growth_transf$P))  #phosphorus concentration
growth_transf$N = scale(log(growth_transf$N))  #nitrogen concentration
growth_transf$la = scale(log(growth_transf$la))  #leaf_area
growth_transf$h = scale(log(growth_transf$h))  #max_height
growth_transf$seed = scale(log(growth_transf$seed))  #seedmass
growth_transf$RGR = scale(log(growth_transf$RGR)) #rgr
growth_transf$vla = scale(log(growth_transf$vla)) #vein length per unit area
growth_transf$wsg = scale(growth_transf$wsg) #wood density
growth_transf$TNCI = scale(log(growth_transf$TNCI)) #total neighbor crowding index
growth_transf$C = scale(growth_transf$C)  #carbon concentration
growth_transf$c13 = scale((growth_transf$c13*-1)) #carbon isotope composition
growth_transf$log.lma = scale(log(growth_transf$lma)) #leaf mass per area

growth_transf = growth_transf[complete.cases(growth_transf),]
rownames(growth_transf) = 1:dim(growth_transf)[1]
 
### M O D E L S  Q #2
m1 = lmer(data = growth_transf, RGR ~ dbh + (1|quad) + (1|species)) #dbh
m2 = lmer(data = growth_transf, RGR ~ log.lma + (1|quad) + (1|species)) #lma
m3 = lmer(data = growth_transf, RGR ~ log_ca2 + (1|quad) + (1|species)) #ca
m4 = lmer(data = growth_transf, RGR ~ mp + (1|quad) + (1|species)) #Mp


Weights(AIC(m1, m2, m3,m5,m6))
c(AIC(m1),AIC(m2),AIC(m3),AIC(m4))

#### K-FOLD CROSS VALIDATION repeated for each model
k = 10 #number of folds
data_cv = growth_transf[sample(nrow(growth_transf)),] #this suffles the data whitout changing the columns
nr = nrow(data_cv)
list_split = split(data_cv, 1:k) #divides the data in k-folds
loss_save_i = c()
###
## model 1
loss.er.Ttry1 = rep(0,k) #saves the loss
for (i in 1:k){
  
  Test = list_split[[i]]
  remove = as.numeric(rownames(Test))
  Train = growth_transf[-remove,]
  
  Ttry1 = m1
  SpredT1 = predict(Ttry1,Test,allow.new.levels = TRUE)
  residual = Test$RGR - SpredT1
  div = round(length(Test$RGR)/10)
  sort = abs(sort(residual,decreasing = T))
  delta = sort[div]
  loss.er.Ttry1[i] = sum(Huber(residual, c = delta))/length(Test$RGR)
  
  rm(Test)
  rm(remove)
  rm(Train)
}
loss_save_i[1] = mean(loss.er.Ttry1)


###########
### For Models using Mpfit all model combinations
##########

g.model = lmer(data = growth_transf, RGR ~ mp+TNCI+vla+c13+P+N+C+wsg+la+h+seed+(1|quad)+(1|species), na.action = "na.fail")

#cluster running
clusterType = if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
cl = try(makeCluster(getOption("cl.cores", 3), type = clusterType))

clusterEvalQ(cl, {
  library(lmerTest)
  library(lme4)
  library(MuMIn)
})
clusterExport(cl, "growth_transf")

#pdredge
drgr1 = pdredge(g.model1,evaluate = T, rank = "AIC",cluster = cl,trace = 2,extra= c(max.r),fixed = "mp2")

## pick the models in which the variables were not highly correlated (multicollinearity <0.6)
max = subset(drgr, max.r < 0.6)
use = (max[1:50,])
cinc = data.frame(use[1:50,])

## models picked from the g.model
m1 = lmer(data = growth_transf, RGR ~ mp + TNCI + wsg +(1|quad) + (1|species))
m2 = lmer(data = growth_transf, RGR ~ mp + TNCI + (1|quad) + (1|species))

w = Weights(c(AIC(m1),AIC(m2)))


#### Prediction graph: predicted vs observed
rgr = unscale(growth_transf$RGR,scale.RGR)
rgr = exp(rgr)

pred.1 = predict(m1, growth_transf, type = "response")
ps1 = unscale(pred.1,scale.RGR)
ps1 = exp(ps1)

plot(ps1,rgr, cex = 2, col = adjustcolor(9,alpha.f=0.05),cex.lab = 1.5,cex.axis=1.5, pch = 16, ylab = "Relative growth rate + 1" , xlab = "Predicted Relative growth rate + 1")
abline(lm(rgr~ps1), col = c("blue"))
cor.test(ps1, rgr,method = "pearson")
legend(0.988,1.24, inset = c(0, 1), c("r = 0.47","p-value <0.001"),bty = "n", xpd = T, plot = T,cex = 1)

## Model plotting

  ci = confint(m2,parm="beta_",method="Wald")
  coefs = data.frame(coef(summary(m1)))
  
  estimate<-c(coefs[3,1],coefs[2,1], coefs[4,1])
  up<-c(coefs[3,1]+1.96*coefs[3,2],coefs[2,1]+1.96*coefs[2,2],coefs[4,1]+1.96*coefs[4,2])
  lo<-c(coefs[3,1]-1.96*coefs[3,2],coefs[2,1]-1.96*coefs[2,2],coefs[4,1]-1.96*coefs[4,2])
  
  l = 3
  odds = estimate 
  index =  l:1
  parameters = c("WD","Mp" ,"NCI") 

  plot(y = index, x = c(-0.5,0,0.5) , cex=0.9, xlab="", ylab=" ", yaxt="n", xaxt="n", xlim=c(-0.5,0.5), ylim=c(0.7, 3.5), type="n" )
  axis(side=2, at=l:1, labels=parameters, cex.axis=0.9)
  axis(side=1, at=c(-0.5,0,0.5), labels=c(-0.5,0,0.5), mgp=c(3,0.3,0), cex.axis=0.9)
  segments(x0=lo, y0=index, x1=up, y1=index)
  abline( v=0, lty=3)
  #abline(h = 1:3,lty=3, lwd=0.5)
  points(y = index, x= odds, cex=0.8, pch=21, bg="white") 
  points(x=odds[c(1,2,3)], y=index[c(1,2,3)], pch=19, cex=0.5, title(xlab="Standarized regression coefficients", line=1.5, cex.lab=0.9))
  
  
 ### function for max multicollinearity 
  max.r <- function(x){
    if(class(x)[length(class(x))] == "lm"){
      corm <- summary(x, correlation=TRUE)$correlation}
    else if(class(x) =="lmerMod"){
      corm <- cov2cor(vcov(x))}
    else if(class(x) =="lmerModLmerTest"){
      corm <- cov2cor(vcov(x))}
    else if(class(x) =="glmerMod"){
      corm <- cov2cor(vcov(x))}
    else if(class(x)=="gls"){
      corm <- summary(x)$corBeta} 
    else if(class(x)=="lme"){
      corm <- summary(x)$corFixed}
    else { print("Error: Invalid model class")}
    corm <- as.matrix(corm)
    if (length(corm)==1){
      corm <- 0
      max(abs(corm))
    } else if (length(corm)==4){
      cormf <- corm[2:nrow(corm),2:ncol(corm)]
      cormf <- 0
      max(abs(cormf))
    } else {
      cormf <- corm[2:nrow(corm),2:ncol(corm)]
      diag(cormf) <- 0
      max(abs(cormf))
    }
  }
  
  