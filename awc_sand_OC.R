mois<-read.csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/Soil_Moisture.csv")
att<-read.csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/Atterberg_limits.csv")
psa<-read.csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/AfSIS_baseline_reference_data/AfSIS_LDPSA.csv")
cn<-read.csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/AfSIS_baseline_reference_data/AfSIS_CN.csv")

mps<-merge(mois,psa)
mpsc<-merge(mps,cn)[,c("SSN","awc1","awc2","psa.c4sand","psa.c4clay","psa.c4silt","Total_Nitrogen","Total_Carbon","Acidified_Carbon")]
str(mpsc)
mpsc<-na.omit(mpsc)
colnames(mpsc)<-c("SSN","awc1","awc2","sand","clay","silt","N","TC","OC")
mpsc$class<-ifelse(mpsc$clay<median(mpsc$clay),"low","high")

m<-with(mpsc,which(awc2==max(awc2)))
n<-with(mpsc,which(awc2==min(awc2)))
s<-with(mpsc,which(sand==max(na.omit(sand))))

mpsc<-mpsc[-c(m,n,s),]
with(mpsc,plot(awc1,OC,xlim=c(0,1)))
with(mpsc,plot(awc1,OC))
with(mpsc,plot(awc2,awc1,xlim=c(0,1)))

pl <-cloud(sand~awc2*OC,data=mpsc)
pl


mpsc<-merge(mps,cn)[,c("SSN","awc1","awc2","psa.c4sand","psa.c4clay","psa.c4silt","Total_Nitrogen","Total_Carbon","Acidified_Carbon")]
str(mpsc)
colnames(mpsc)<-c("SSN","awc1","awc2","sand","clay","silt","N","TC","OC")

m<-with(mpsc,which(awc2==max(awc2)))
n<-with(mpsc,which(awc2==min(awc2)))
s<-with(mpsc,which(clay==max(na.omit(clay))))

mpsc<-mpsc[-c(m,n,s),]
with(mpsc,plot(awc1,OC,xlim=c(0,1)))
with(mpsc,plot(awc1,OC))
with(mpsc,plot(awc2,awc1,xlim=c(0,1)))

pl <-cloud(awc2~clay*OC,data=mpsc)
pl