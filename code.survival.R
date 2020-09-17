@@ -0,0 +1,248 @@
### SURVIVAL ####
##############
rm(list = ls())
set.seed(NULL)
packages = c("lmerTest","lme4","MuMIn","qrmix","data.table","MLmetrics","parallel")
lapply(packages,require,character.only=T)

##---
survival = fread("~/Data.txt")
initial = as.data.frame(survival)
survival_transf = subset(initial, !is.na(ca2)) #canopy
survival_transf = subset(survival_transf,!is.na(survival))
survival_transf = subset(survival_transf,TNCI >0 ) #total neighborhood index

survival_transf$RGR = NULL
survival_transf$log_ca2 = log(survival_transf$ca2)
survival_transf$log_ca2 = scale(survival_transf$log_ca2) 

survival_transf$mp = log(survival_transf$ca2*survival_transf$lma) 
survival_transf$mp = scale(survival_transf$mp)
survival_transf = survival_transf[complete.cases(survival_transf),]

#Log transformation and/or scale
survival_transf[,10:19] = scale(log(survival_transf[,10:19]))
survival_transf$C = scale(survival_transf$C)  
survival_transf$c13 = scale(log(survival_transf$c13*-1))
survival_transf$wsg = scale(survival_transf$wsg)
rownames(survival_transf) = 1:dim(survival_transf)[1]


### M O D E L S  Q #2
sm1 = glmer(data = survival_transf, survival ~ dbh + (1|quad) + (1|species),family=binomial(link=logit)) #dbh
sm2 = glmer(data = survival_transf, survival ~ lma + (1|quad) + (1|species),family=binomial(link=logit)) #lma
sm3 = glmer(data = survival_transf, survival ~ log_ca2 + (1|quad) + (1|species),family=binomial(link=logit)) #ca
sm4 = glmer(data = survival_transf, survival ~ mp+ (1|quad) + (1|species),family=binomial(link=logit)) #Mp

Weights(AIC(sm1, sm2, sm3,sm4))
c(AIC(sm1),AIC(sm2),AIC(sm3),AIC(sm4))

#### K-FOLD CROSS VALIDATION repeated for each model
set.seed(25)
k = 10 #number of folds
data_cv = survival_transf[sample(nrow(survival_transf)),] #this suffles the data whitout changing the columns
nr = nrow(data_cv)
list_split = split(data_cv, 1:k) #divides the data in k-folds
loss_save_i = c()

## model 1
loss.er.Ttry1 = rep(0,k) #saves the loss
for (i in 1:k){
  
  Test = list_split[[i]]
  remove = as.numeric(rownames(Test))
  Train = survival_transf[-remove,]
  
  Ttry1= sm1
  SpredT1 = predict(Ttry1,Test,type = "response")
  loss.er.Ttry1[i] = LogLoss(SpredT1,Test$survival)
  
  rm(Test)
  rm(remove)
  rm(Train)
}
loss_save_i[1] = mean(loss.er.Ttry1)

###########
## fit all model combinations
### For Models using Mp (mp) and CA (log_ca2)

global.model = glmer(data = survival_transf, survival~mp+vla+c13+P+N+C+wsg+la+h+seed+(1|quad)+(1|species),family=binomial(link=logit), glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
,na.action = "na.fail")

global.model2 = glmer(data = survival_transf, survival~log_ca2+TNCI+vla+c13+P+N+C+wsg+la+h+seed+(1|quad)+(1|species),family=binomial(link=logit), glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
,na.action = "na.fail")

#cluster running
clusterType = if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
cl = try(makeCluster(getOption("cl.cores", 3), type = clusterType))

clusterEvalQ(cl, {
  library(lmerTest)
  library(lme4)
  library(MuMIn)
})
clusterExport(cl, "survival_transf")

#pdredge for mp
dr.s.mp = pdredge(global.model, cluster = cl, beta = c("none"),evaluate = T, rank = "AIC", trace = 2, extra= c(max.r), fixed = "mp")
dr.s.ca = pdredge(global.model2, cluster = cl, beta = c("none"),evaluate = T, rank = "AIC", trace = 2, extra= c(max.r), fixed = "log_ca2")

## pick the models in which the variables were not highly correlated (multicollinearity <0.6)
mp.max = subset(dr.s.mp, max.r< 0.6)
mp.use = (mp.max[1:52,])
mp.cinc = data.frame(mp.use[1:52,])

ca.max = subset(dr.s.ca, max.r< 0.6)
ca.use = (ca.max[1:50,])
ca.cinc = data.frame(ca.use[1:50,])

mp.delta = data.frame(mp.max[1:21,])
ca.delta = data.frame(ca.max[1:7,])

## models picked from the global.model using Mp
am1 = glmer(data = survival_transf, survival~mp+C+h+P+(1|quad)+(1|species),family=binomial(link=logit))
am2 = glmer(data = survival_transf, survival~mp+C+h+P+N+seed+vla+(1|quad)+(1|species),family=binomial(link=logit)) 
am3 = glmer(data = survival_transf, survival~mp+C+h+P+TNCI+(1|quad)+(1|species),family=binomial(link=logit)) 
am4 = glmer(data = survival_transf, survival~mp+C+h+P+N+seed+TNCI+vla+(1|quad)+(1|species),family=binomial(link=logit)) 
am5 = glmer(data = survival_transf, survival~mp+P+(1|quad)+(1|species),family=binomial(link=logit))
am6 = glmer(data = survival_transf, survival~mp+C+P+(1|quad)+(1|species),family=binomial(link=logit))
am7 = glmer(data = survival_transf, survival~mp+C+h+P+c13+(1|quad)+(1|species),family=binomial(link=logit))
am8 = glmer(data = survival_transf, survival~mp+P+TNCI+(1|quad)+(1|species),family=binomial(link=logit)) 
am9 = glmer(data = survival_transf, survival~mp+C+h+P+N+seed+vla+c13+(1|quad)+(1|species),family=binomial(link=logit)) 
am10 = glmer(data = survival_transf, survival~mp+C+P+TNCI+(1|quad)+(1|species),family=binomial(link=logit)) 
am11 = glmer(data = survival_transf, survival~mp+h+P+(1|quad)+(1|species),family=binomial(link=logit))
am12 = glmer(data = survival_transf, survival~mp+C+P+wsg+(1|quad)+(1|species),family=binomial(link=logit))
am13 = glmer(data = survival_transf, survival~mp+C+h+P+TNCI+c13+(1|quad)+(1|species),family=binomial(link=logit))
am14 = glmer(data = survival_transf, survival~mp+C+h+P+wsg+(1|quad)+(1|species),family=binomial(link=logit)) 
am15 = glmer(data = survival_transf, survival~mp+h+P+TNCI+(1|quad)+(1|species),family=binomial(link=logit)) 
am16 = glmer(data = survival_transf, survival~mp+P+la+(1|quad)+(1|species),family=binomial(link=logit))
am17 = glmer(data = survival_transf, survival~mp+C+P+TNCI+wsg+(1|quad)+(1|species),family=binomial(link=logit))
am18 = glmer(data = survival_transf, survival~mp+C+c13+(1|quad)+(1|species),family=binomial(link=logit)) 

## models picked from the global.model2 using CA
m1 = glmer(data = survival_transf, survival~log_ca2+C+h+P+(1|quad)+(1|species),family=binomial(link=logit))
m2 = glmer(data = survival_transf, survival~log_ca2+C+h+P+TNCI+(1|quad)+(1|species),family=binomial(link=logit)) 
m3 = glmer(data = survival_transf, survival~log_ca2+C+h+N+P+seed+vla+(1|quad)+(1|species),family=binomial(link=logit)) 
m4 = glmer(data = survival_transf, survival~log_ca2+C+h+N+P+seed+TNCI+vla+(1|quad)+(1|species),family=binomial(link=logit)) 
m5 = glmer(data = survival_transf, survival~log_ca2+C+h+P+c13+(1|quad)+(1|species),family=binomial(link=logit)) 
m6 = glmer(data = survival_transf, survival~log_ca2+C+h+P+wsg+(1|quad)+(1|species),family=binomial(link=logit)) 


##############################################################################################################
##############################################################################################################
#### Model Averaging
##############################################################################################################
#######
#### Mp

Weights(AIC(am1, am2, am3, am4,am5,am6,am7, am8, am9, am10, am11, am12,am13, am14, am15, am16, am17, am18)) 
lista1 = list(am1,am2,am3,am4,am5,am6,am7, am8, am9, am10, am11, am12,am13, am14, am15, am16, am17, am18) 
av1 = model.avg(lista1)
summary(av1)

##############################################################################################################
#######
#### CA
Weights(AIC(m1, m2, m3, m4)) #,m15,m16))
lista2 = list(m1, m2, m3, m4,m5,m6) #,am15,am16,am17,am18,am19,am20)
av2 = model.avg(lista2)
summary(av2)

# average model plotting Mp

par(mfrow = c(1,2))
estimates = c(0.862657, -0.426770,-0.414761,0.298176,0.086368,-0.073484,0.053256,-0.035824,0.035221,0.010233,0.009505)
se = c(0.028991446,0.354798375,0.244282467,0.348005285,0.159057061,0.135574247,0.109038205,0.079246190,0.073046289,0.016816454,0.022248649)
up = estimates+(1.96*se)
lo = estimates-(1.96*se)
l = 11
index =  l:1
parameters = c("Mp","C", "P" , "Hmax" , "N","SM" , "VLA", "WD", "c13", "NCI","LA") 

### plot
plot(y = 5:1, x = c(-1.3,-0.5,0,0.5,1.3) , cex=1, xlab="", ylab=" ", yaxt="n", xaxt="n", xlim=c(-1.4,1.4), ylim=c(0.7, 11.5), type="n" )
axis(side=2, at=l:1, labels=parameters, cex.axis=1)
axis(side=1, at=c(-1,-0.5,0,0.5,1), labels=c(-1,-0.5,0,0.5,1), mgp=c(3,0.3,0), cex.axis=1)
segments(x0=lo, y0=index, x1=up, y1=index)
abline( v=0, lty=3)
#abline(h = 1:3,lty=3, lwd=0.5)
points(y = index, x= estimates, cex=1, pch=21, bg="white") 
points(x=estimates[c(1)], y=index[c(1)], pch=19, cex=1, title(xlab="Standarized regression coefficients", line=1.5, cex.lab=1))

# model plot CA
estimates = c(0.794230, -0.649450,0.527830,-0.500030,0.119970,-0.102100,0.101960,0.012350,0.010424,-0.008646)
se = c(0.026652547,0.295509256,0.247347906,0.206352467,0.183109650,0.156147608,0.156818168,0.017786665,0.030456859,0.036400774)
up = estimates+(1.96*se)
lo = estimates-(1.96*se)
l = 10
index =  l:1
parameters = c("CA","C" , "Hmax" ,"P", "N" , "SM","VLA", "NCI","c13","WD") 

### plot
plot(y = 5:1, x = c(-1.3,-0.5,0,0.5,1.3) , cex=1, xlab="", ylab=" ", yaxt="n", xaxt="n", xlim=c(-1.4,1.4), ylim=c(0.7, 10.5), type="n" )
axis(side=2, at=l:1, labels=parameters, cex.axis=1)
axis(side=1, at=c(-1,-0.5,0,0.5,1), labels=c(-1,-0.5,0,0.5,1), mgp=c(3,0.3,0), cex.axis=1)
segments(x0=lo, y0=index, x1=up, y1=index)
abline( v=0, lty=3)
#abline(h = 1:3,lty=3, lwd=0.5)
points(y = index, x= estimates, cex=1, pch=21, bg="white") 
points(x=estimates[c(1,2,3,4)], y=index[c(1,2,3,4)], pch=19, cex=1, title(xlab="Standarized regression coefficients", line=1.5, cex.lab=1))

#### Prediction graph: predicted vs observed
library(verification)
aveg.1 = predict(av1,survival_transf, type = "response",full=T)
survival_transf$ave1 = aveg.1
roc.area(survival_transf$survival,survival_transf$ave1)$A

aveg.2 = predict(av2,survival_transf, type = "response",full=T)
survival_transf$ave2 = aveg.2
roc.area(survival_transf$survival,survival_transf$ave2)$A

par(mfrow = c(1,2), mai = c(1,1.5,1,1)) #bottom, left, top, right)
roc.plot(survival_transf$survival,aveg.1, binormal = TRUE, plot.thres = NULL, main = " Mp - ROC Curve
AUC = 0.765", ylab = "Sensitivity", xlab = "1 - Specificity", lwd = 4,cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)

roc.plot(survival_transf$survival,aveg.2, binormal = TRUE, plot.thres = NULL, main = "CA - ROC Curve
AUC = 0.764", ylab = "Sensitivity", xlab = "1 - Specificity", lwd = 4,cex=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)


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

