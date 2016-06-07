suppressMessages(library(caret))
suppressMessages(library(readr))
suppressMessages(library(mlbench))
suppressMessages(library(pROC))
suppressMessages(library(rpart))
suppressMessages(library(caretEnsemble))
suppressMessages(library(soil.spec))
suppressMessages(library(FitAR))
suppressMessages(library(doParallel))
suppressMessages(library(ggplot2))
registerDoParallel()
getDoParWorkers()
require(readr)
#boxcox transformation for response variables
#1.MSC
#2.SNV
#3.snv-detrend
#Read notes 1 and 2 for correct computation of r.squared values.
a<-"~/Models"
setwd(a)
#read first derivative spectra
suppressMessages(library(soil.spec))

mir<-read_csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/Calibration_Htsxt_MIR.csv")
ref<-read_csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/AfSIS_baseline_reference_data/AfSIS_soil_moisture.csv")#[,-c("SSN","awc1","awc2")]
ref<-na.omit(ref[,c("SSN","awc1","awc2")])

m<-which(ref[,"awc1"]==max(na.omit(ref[,"awc1"]))) #Extremely high
n<-which(ref[,"awc1"]==min(na.omit(ref[,"awc1"]))) #Extremely low
ref<-ref[-c(m,n),]#Remove low and extreme awc1 and awc2 values.
#preprocess by first derivative
#Set the spectra table as a matrix
mim<-as.matrix(mir[,-1])

#Set the colnames of mim matrix as numeric
colnames(mim)<-as.numeric(substr(colnames(mim),2,15))

#Preprocess the raw spectra by first derivative; use library soil.spec
de<-trans(mim,tr="derivative",order=1,gap=21)
der<-as.matrix(de$trans)

dc<-ncol(der)

#Combine derivative with the ssn
derssn<-cbind(as.vector(mir[,1]),der)
colnames(derssn)<-c("SSN",colnames(der))
write.table(derssn, file="AfSIS_first_derivative_mir_data.csv",sep=",",row.names=FALSE)
mir<-read_csv("AfSIS_first_derivative_mir_data.csv")
#Merge preprocessed data with reference
cal.0<-merge(ref,mir,by="SSN")

##### Here check the variable names to know where reference starts and ends and the same for IR.
hd<-colnames(cal.0)[2:3]# 6:7 are the columns with maize yld data; this will change for different data set
b<-getwd()
if(!file.exists("Models")){dir.create("Models")}
if(!file.exists("calibration_plots")){dir.create("calibration_plots")}

pms.a <-NULL
pred.all <-NULL
for (q in 1:length(hd)){
cal<-cal.0[,c(hd[q],colnames(mir[,-1]))]
set.seed(1009)
#Select training and testing sets
cal<-na.omit(cal)

trainX <-cal[, -1]
trainY <-cal[,1]

#train <- createFolds (trainY, k=round(0.2*length(trainY),)) Not efficient picks fw smples every time
#train <-createMultiFolds(trainY,k=10,time=3) time consuming
train <-createDataPartition(trainY,time=3,p=0.5)
#mir.grid<-expand.grid(degree =1, nprun = (1:20)*3)


plsrf <- train(trainY~., method="rf", data=cal[,-1],
			tuneLength =3,
			metric="RMSE",
			maximise=FALSE,
			.mtry=50,
			trControl=trainControl(method="oob", indexOut = train,p=0.75))
			#,tuneGrid=mir.grid)
		
trees<-as.numeric(plsrf$bestTune)[1]
l<- which(plsrf$results$mtry==trees)

#Note1: Get final model to compute coefficient for variation explained
predi<-plsrf$finalModel$predicted
y<-plsrf$finalModel$y
rsq<-round(1-sum((y-predi)^2)/sum((y-mean(y))^2),2)
pred<-predict(plsrf,trainX)
#Note2:this correlation coefficient is variance explained
#cor(pred,trainY)^2
#rsq<-round(cor(pred,trainY)^2,3)
round(summary(lm(pred~trainY))$r.squared,3)
rmse<-round(plsrf$results$RMSE[l],2)
setwd(paste(b,"Models",sep="/"))#Save fitted model
#Get model summary
pms<-c(hd[q],rsq,rmse,trees)
pms.a <- rbind(pms.a,pms)
saveRDS(plsrf,file=paste(hd[q],".Rds",sep=""))

#Get predicted values then take inverse bxcx
colnames(mir)<-c("SSN",names(plsrf$trainingData))
pred.a<-predict(plsrf,mir[,-1])
pred.all<-cbind(pred.all,as.vector(pred.a))
setwd(paste(b,"calibration_plots",sep="/"))#Save PLS models
png(file=paste0(hd[q],".png"))
plot(pred,trainY,ylim=c(min(pred,trainY),max(pred,trainY)),xlim=c(min(pred,trainY),max(pred,trainY)),ylab="",xlab="",pch=19,col="blue")
abline(a=0,b=1,col="red")
require(stats)
#r<-round(cor(exp(pred$fitted.values[,,nc]),na.omit(cutchem[,q]),method="pearson")^2,2)
mtext(paste(hd[q]," ","n=",nrow(cal),sep=""),side=3,line=0.7,cex=1.8)
mtext(paste("Predicted ",sep=""),side=1,line=2.4,cex=1.6)
mtext(paste("Measured ",sep=""),side=2,line=2.4,cex=1.6)
mylabel = bquote(bold(r)^2==.(format(rsq,digits=3)))
#legend("topleft",c(paste("r.squared",rsq,sep="="),paste("rmse=",rmse,sep="")),bty="n")
ab<-mean(c(pred,trainY))
bb<-max(pred)
text(x=ab,y=bb,labels =mylabel,cex=1.3)
dev.off()
}
#Save pms
colnames(pms.a)<- c("Property","R_sq.","RMSE","Mtry")
write.table(pms.a,file="Model_Summary.csv",sep=",",row.names=FALSE)

#Save predicted values
#Add SSN
pred.all<-cbind(as.vector(mir[,1]),pred.all)
colnames(pred.all) <- c("SSN", hd)
setwd(b)
write.table(pred.all,file="All Predictions.csv",sep=",",row.names=FALSE)

