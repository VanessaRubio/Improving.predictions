###### Allometric graphs

data = species.allom
data$dbh = log(data$dbh,base = 10)
data$AVE.RADIUS = log(data$AVE.RADIUS,base = 10)
#data$area = pi*((data$AVE.RADIUS)^2)
#data$area = log10(data$area)

spp = sort(unique(data$spp))
results = data.frame(spp,rep(0,17),rep(0,17),rep(0,17),rep(0,17),rep(0,17),rep(0,17))
colnames(results) = c("spp","n","intercept","slope","r2","st.error.int","st.error.slope")

par(mfrow = c(5,4),mai = c(0.62,0.62,0.4,0.2))
for (i in 1:17){
  
  sp = spp[i]
  tmp = data[data$spp==sp,]
  results[i,"n"] = nrow(tmp)
  l.m = lm(tmp$AVE.RADIUS~tmp$dbh)
  plot(tmp$dbh,tmp$AVE.RADIUS,cex=3.5,xlab = "log(dbh)",ylab = "log(crown radius)",ylim = c(-0.8,1),xlim = c(-0.4,1.8),cex.axis = 1.3,cex.lab = 1.7, col = adjustcolor(9,alpha.f=0.5), pch = 20)
  abline(lm(tmp$AVE.RADIUS~tmp$dbh))
  title(paste(sp,", n =",nrow(tmp)),cex.main = 1.5)
  results[i,"intercept"] = l.m$coefficients[1]
  results[i,"slope"] = l.m$coefficients[2]
  results[i,"st.error.int"] = (summary(l.m))$coefficients[3]
  results[i,"st.error.slope"] = (summary(l.m))$coefficients[4]
  results[i,"r2"] = (summary(l.m))$adj.r.squared
}
