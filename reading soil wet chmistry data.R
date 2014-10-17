soil<-read.csv(file.choose())
soil<-wharevername
dim(soil)

attach(soil)

colnames(soil)

ph<-mean(PH) 
round(ph,1)

phg<-subset(soil,soil$PH>4.8)
phg

setwd("e:/R")
getwd()#To know your current wd
write.table(phg,file="Weak acidic soil.csv",sep=",",row.names=F)
