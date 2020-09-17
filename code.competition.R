rm(list = ls())

library(data.table)

c05 = read.table("~census2005.txt",h=T)
colnames(c05) = c("species","quad","lxm","lym","TreeID","Tag","StemID","StemTag","dbh")
c05$TNCI = rep(0,dim(c05)[1]) #Total Neighbor crowding index
c05$dbh = c05$dbh/100 #transform to meters

### Total Neighbor crowdind index (TNCI)

calculate_tnci = function(x){
  
  r = 20 #radius for analysis 20m
  use_row = as.numeric(rownames(x))
  others = c05[-use_row,]  #finds the other trees removing the row with the focal
  rownames(others) = 1:dim(others)[1] #changes the rownames so the positions are correct
  distances = sqrt((others[,"lxm"]-x[,"lxm"])^2 + (others[,"lym"]-x[,"lym"])^2 ) #calculated distances
  positions = which(distances<=r) #less than 20m of radius
  #all_sp = as.character(others[positions,"species"])  #species names of the trees in the radius
  #tags = others[positions,"Tag"]  #saves the tags ids of the trees in a 20m radius
  dbhs = others[positions,"dbh"]  #saves the dbhs of the trees in a 20m radius
  ds = distances[positions] #distances from each tree to the focal tree
  ds[which(ds==0)] = 0.3 #if the three was positioned in the census in the same xy coordintes, I assumed 30cm in distance, otherwise Inf and removes the tree

  tnci = sum((dbhs)^2/(ds)^2,na.rm = T)
  
  return(tnci)
}

for (i in 1:dim(c05)[1]){
  use = c05[i,] #focal tree
  c05[i,"TNCI"] = calculate_tnci(use)
  print(i)
}

c05$dbh = c05$dbh*100 #restarts the initial dbh for these analysis to mm

