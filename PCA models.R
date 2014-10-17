library(soil.spec)

#Read Ghana -AfSIS samples
#Read AfSIS spectra
afs<-read.csv("~/GRP4/Studies_data/Job_Kihara/Crop_response/Data/Averaged DT MIR.csv")
yld<-read.csv("~/GRP4/Studies_data/Job_Kihara/Crop_response/Data/Crop_response.csv")
str(yld)
afsn<-afs[,-c(1,3579)]
colnames(afsn)<-as.numeric(substr(colnames(afsn),2,16))

#Obtain derivatives
de.s<-trans(afsn,tr="derivative",order=1,gap=21)
der.s<-as.matrix(de.s$trans)


#########Run PCA#################################
pc<-prcomp(der.s)
summary(pc)
imp<-summary(pc)$importance
pcs<-pc$x[,1:10]
pcs.ssn<-cbind(as.vector(afs[,1]),pcs)
colnames(pcs.ssn)<-c("SSN",colnames(pcs))
ypc<-merge(yld,pcs.ssn)
write.table(ypc,file="~/GRP4/Studies_data/Job_Kihara/Crop_response/Data/pc scores plus codes.csv",sep=",",row.names=FALSE)
ypc<-read.csv("~/GRP4/Studies_data/Job_Kihara/Crop_response/Data/pc scores plus codes.csv")
ypc.c<-subset(ypc,ypc$TrtDesc=="Control")
ypc.cs<-ypc.c[,-c(1:ncol(yld))]
ypc.n<-subset (ypc,ypc$TrtDesc!="Control")
ypc.ns<-ypc.n[,-c(1:ncol(yld))]
ypc.cs<-as.data.frame(ypc.cs)
png(file="~/GRP4/Studies_data/Job_Kihara/Crop_response/Soil spectral PCA.png",width=600,height=600)

plot(ypc.cs[,1:2],col="green4",pch=19,main="PCA scores for soil spectra",xlab=paste0("PC1 explains ", round(imp[2,1],3)*100, " % total variation"),ylab=paste0("PC1 explains ", round(imp[2,2],3)*100, " % total variation"))
points(ypc.ns[,1:2],col="yellow2",pch=19)
legend("bottomright",pch=19,col=c("green4","yellow2"),c("Control","Treated"),bty="n",text.col=c("green4","yellow2"))
dev.off()
