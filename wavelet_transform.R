library(wmtsa)
library(wavelets)
library(soil.spec)
library(ChemoSpec)
ir<-read.csv(file.choose())
#mir<-ir[,-c(2:1813)]
setwd("~/Thesis/Data/wavelets")
str(mir)
par(mfrow=c(2,2))
for ( i in 1:4){
	ir.i.all<-NULL
for(k in 1:nrow(mir)){
t2<-wavDWT(mir[k,-c(1,ncol(mir))],wave="haar",keep.series=FALSE)
names(t2)
par(mfrow=c(1,1))
j<-length(t2$data)-2
minnn<-NULL
 for ( a in 1:j){
 	minn<-min(t2$data[[a]])
 	minnn<-c(minnn,minn)
 	}
 maxxx<-NULL	
 for ( b in 1:j){
 	maxx<-max(t2$data[[b]])
 	maxxx<-c(maxxx,maxx)
 	}

plot(t2$data[[i]],col=i,type="l")
ir.i<-t2$data[[i]]
ir.i.all<-rbind(ir.i,ir.i.all)}
ir.i.all<-cbind(as.vector(ir[,1]),ir.i.all)
colnames(ir.i.all)<-c("SSN",paste0("D",i,".",2:ncol(ir.i.all)-1))
write.table(ir.i.all,file=paste0("Wavelet transform level_",i,".csv"),sep=",",row.names=FALSE)	

}
plot(t2$data$d1,type="l",ylim=range(minnn,maxxx))
for ( i in 2:j){
lines(t2$data[[i]],col=i)}