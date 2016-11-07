###  Script for doing calibration models using IR data
# Andrew Sila , April 2016
#Begin by checking the required packages are installed
#####
#is.installed<-function(anypkg){
	#is.element(anypkg, installed.packages()[,1])
#}

#if (!is.installed(c("hexView1","hexView11","hexView")){
#paste("The following packages will be installed: ", c("hexView1","hexView11","hexView")[which(!is.installed(c("hexView1","hexView11","hexView"))==TRUE)])
#}

#install.answer<-function(){
	#n<-readline(prompt="Install the packages?:")
	#return()
#}
#if (!is.installed("hexView")){
	#install.packages("hexView")
#}	
### Load packages
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
suppressMessages(library(dplyr))
suppressMessages(library(downloader))
suppressMessages(library(prospectr))
suppressMessages(library(randomForest))
suppressMessages(library(gridExtra))

registerDoParallel()
getDoParWorkers()

### To install any required pacakge use:
#install.packages(c("package1, package2,...packagen"), dependencies=TRUE)

#Start by setting working directory
calibrate<-function(wd,infrared.data,reference.data,hout,method=c("RF","PLS")){
	if(method=="RF"){
#a<-wd#Mac OSX
#a<-"D:/oaf"#Windows
setwd(wd)
#download("https://www.dropbox.com/s/em8w2cdnnnkdrvo/rwanda.zip","rwanda.zip", mode="wb")
#unzip("rwanda.zip",overwrite=T)

#read raw IR data then preprocess by first derivatives using SG algorithm.
#mir<-read_csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/Calibration_Htsxt_MIR.csv")
mir<-infrared.data
ref<-reference.data
hout<-hout
#Read reference data
#ref<-read_csv("~/Studies_data/data/Calibration_ref_data.csv")
#k<-which(mir$SSN%in%ref$ICRAF_ID)

#mir<-mir[k,]
#Step 2: Preprocess with SG algorithm
mir1<-as.matrix(mir[,-1])#Exclude metadata variables
wave<-as.numeric(substr(colnames(mir1),2,19))
colnames(mir1)<-wave
#First derivative
de1<-trans(mir1,tr="derivative",order=1,gap=23)
der1<-rev(as.data.frame(de1$trans))
colnames(der1)<-paste0("m",wave)

#Save derivative spectra
der1.ssn<-as.data.frame(cbind(as.vector(mir[,1]),der1)); colnames(der1.ssn)<-c("SSN",colnames(der1))
write.table(der1.ssn,file="First derivative.csv",sep=",",row.names=FALSE)
der1.ssn<-as.data.frame(read_csv("First derivative.csv"))
#merge with first derivative preprocessed spectra
ref.mir<-merge(ref,der1.ssn,by.x="SSN",by.y="SSN")
rc<-colnames(ref)
#which columns contains reference data?
ref<-ref.mir[,rc]
#Extract spectral predictors
mirp<-colnames(der1.ssn)[-1]
spectra<-der1.ssn[,mirp]
#Create two new subfolders within the current working using:
b<-getwd()
if(!file.exists("Models")){dir.create("Models")}
if(!file.exists("calibration_plots")){dir.create("calibration_plots")}

#Fit calibration models for the training set and use the testing set to validate the models.
#which are the training samples?
set.seed(67523)
#testing<-sample(1:nrow(ref.mir),0.33*nrow(ref.mir))#Hold-out a third of the calibration set for validation
#testing<-which(ref.mir$set=="val") #Use the val set rows
testing<-which(der1.ssn$SSN%in%hout$SSN)#with hout

#Use Kennard_Stone
sel <- kenStone(spectra,k=round(0.7*nrow(spectra)),pc=.99)
#View selected samples
plot(sel$pc[,1:2],xlab='PC1',ylab='PC2')
points(sel$pc[sel$model,1:2],pch=19,col=2) # points selected for calibration
#testing<-sel$model

#### This chunk is optional; run it if you want to share the testing and training set with other people################
#Set validation rows to 2 and calibration to 1 and save in a new column -set
#ref.mir$set<-rep(NA,nrow(ref.mir))
#ref.mir$set[testing]<-2
#ref.mir$set[-testing]<-1
#p<-which(colnames(ref.mir)=="set")
#Save the file with the new column set
#Run the next line to round of the wavenumbers into one dp; 
#colnames(ref.mir)<-c(colnames(ref.mir[,1:20]),paste0("m",round(as.numeric(colnames(ref.mir[,21:1782])),1)),"Set")
#write.table(ref.mir,file="~/Studies_data/Rwanda/RF_models/Derivative_Calibration_Holdout.csv",sep=",",row.names=FALSE)
#### End of optional chunk################

#Loop for calibration of all soil properties in the reference set starts here
msummary<-NULL
hd<-colnames(ref)[-1]
for (q in 1:length(hd)){
refq<-which(colnames(ref)%in%hd[q])
ref.q<-ref[,refq]
pms.a <-NULL
pred.all <-NULL
cal<-cbind(as.vector(ref.q),spectra)[-testing,]
val<-cbind(as.vector(ref.q),spectra)[testing,]
colnames(cal)<-c(colnames(ref)[refq],colnames(spectra))
colnames(val)<-colnames(cal)
cal<-na.omit(cal)
val<-na.omit(val)
trainX <-cal[, -1]
trainY <-cal[,1]
#mgrid<-expand.grid(.mtry=seq(50,600,50))
#rf.m<-randomForest(trainY~.,data=cal[,-1],corr.bias=TRUE)
rf.m <- randomForest(data=cal, x=cal[,-1], y=cal[,1], ntree=200, importance=TRUE, na.action=na.omit)


predicted<-predict(rf.m )
measured<-trainY
training.parameters<-round(postResample(predicted,measured),2)#computes RMSE and R-squared values for the calibration set
RSQ<-training.parameters[2]
RMSE<-training.parameters[1]
predicted.test<-predict(rf.m,val[,-1])
measured.test<-val[,1]

testing.parameters<-round(postResample(predicted.test,measured.test),2)#computes RMSE and R-squared values for the validation set
RSP<-testing.parameters[2]
RMSEP<-testing.parameters[1]
model.summary<-c(hd[q],training.parameters,testing.parameters)
msummary<-rbind(msummary,model.summary)
saveRDS(rf.m,file=paste0(b,"/","models/",hd[q],".rds"))

pm<-as.data.frame(cbind(measured,predicted))

p<-ggplot(pm, aes(x=measured,y=predicted))+
geom_point(col="black",size=3,alpha=0.2)+
ggtitle(paste0("Calibration for ",hd[q]))+
	#xlim(as.numeric(as.vector(substr(variable,2,19))))+
	xlab("Measured")+
	ylab("Predicted")
	#theme with white background
  #theme_bw() +

  #eliminates background, gridlines, and chart border
  #theme(
   # plot.background = element_blank()
   #,panel.grid.major = element_blank()
   #,panel.grid.minor = element_blank()
  #)
p<-p+stat_smooth(method=lm, se=FALSE, color='black',alpha=0.1)
p<-p+theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20))
p <- p + theme(text = element_text(size = 20)) # this will change all text size 
#p<-p+theme(axis.title.x = element_text(size = rel(2)))+theme(axis.title.y = element_text(size = rel(2)))
p<-p + annotate('text', label=paste('R^2==',RSQ), parse=TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8) +annotate('text', label=paste('RMSE==',RMSE), parse=TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)

pmp<-as.data.frame(cbind(measured.test,predicted.test))#Validation data
p2<-ggplot(pmp, aes(x=measured.test,y=predicted.test))+
geom_point(col="brown",size=3,alpha=0.2)+
ggtitle(paste0("Validation for ",hd[q]))+
	#xlim(as.numeric(as.vector(substr(variable,2,19))))+
	xlab("Measured")+
	ylab("Predicted")
	#theme with white background
  #theme_bw() +

  #eliminates background, gridlines, and chart border
  #theme(
   # plot.background = element_blank()
   #,panel.grid.major = element_blank()
   #,panel.grid.minor = element_blank()
  #)
p2<-p2+stat_smooth(method=lm, se=FALSE, color='brown',alpha=0.1)
p2<-p2+theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20))
p2 <- p2 + theme(text = element_text(size = 20)) # this will change all text size 
#p<-p+theme(axis.title.x = element_text(size = rel(2)))+theme(axis.title.y = element_text(size = rel(2)))
p2<-p2 + annotate('text', label=paste('R^2==',RSP), parse=TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8) +annotate('text', label=paste('RMSE==',RMSEP), parse=TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)

######################################

png(file=paste0(b,"/Calibration_plots/",hd[q],".png"),height=400,width=800)
grid.arrange(p,p2,nrow=1)
dev.off()
}
colnames(msummary)<-c("Soil properties","OOB RMSEC","OOB Rsquared", "Holdout RMSEP","Holdout Rsquared")
write.table(msummary,file="Model_Summary.csv",sep=",",row.names=FALSE)

### All Samples #####
b<-getwd()
if(!file.exists("Full_Models")){dir.create("Full_Models")}
if(!file.exists("Full_calibration_plots")){dir.create("Full_calibration_plots")}
msummary<-NULL
hd<-colnames(ref[,-1])#Exclude SSN 
all.predicted<-NULL
for (q in 1:length(hd)) {
refq<-which(colnames(ref)%in%hd[q])
ref.q<-ref[,refq]
pms.a <-NULL
pred.all <-NULL
cal<-cbind(as.vector(ref.q),spectra)
colnames(cal)<-c(colnames(ref)[refq],colnames(spectra))
p<-which(is.na(der1.ssn[,1])==TRUE)
ifelse(length(p)>0,ssn<-der1.ssn[-p,1],ssn<-der1.ssn[,1])
ifelse(length(p)>0,der1.ssn<-der1.ssn[-p,],der1.ssn<-der1.ssn)

#Select training and testing sets
cal<-na.omit(cal)
trainX <-cal[, -1]
trainY <-cal[,1]

#rf.m<-randomForest(trainY~.,data=cal[,-1],corr.bias=TRUE,metric=metric, tuneGrid=tunegrid, trControl=control)
rf.m <- randomForest(data=cal, x=cal[,-1], y=cal[,1], ntree=200, importance=TRUE, na.action=na.omit)
predi<-predict(rf.m)
y<-trainY
training.parameters<-c(hd[q],round(postResample(predi,y),3))
msummary<-rbind(msummary,training.parameters)
saveRDS(rf.m,file=paste0(b,"/","Full_Models/",hd[q],".rds"))
#Training
predicted<-predict(rf.m )
measured<-trainY
training.parameters<-round(postResample(predicted,measured),3)#computes RMSE and R-squared values for the calibration set
RSQ<-training.parameters[2]
RMSE<-training.parameters[1]
predicted.test<-predict(rf.m,val[,-1])
measured.test<-val[,1]

testing.parameters<-round(postResample(predicted.test,measured.test),2)#computes RMSE and R-squared values for the validation set
RSP<-testing.parameters[2]
RMSEP<-testing.parameters[1]
model.summary<-c(hd[q],training.parameters,testing.parameters)
msummary<-rbind(msummary,model.summary)
saveRDS(rf.m,file=paste0(b,"/","models/",hd[q],".rds"))

pm<-as.data.frame(cbind(y,pred))
colnames(pm)<-c("measured","predicted")
p<-ggplot(pm, aes(x=measured,y=predicted))+
geom_point(col="black",size=3,alpha=0.2)+
ggtitle(paste0("Calibration for ",hd[q]))+
	#xlim(as.numeric(as.vector(substr(variable,2,19))))+
	xlab("Measured")+
	ylab("Predicted")
	#theme with white background
  #theme_bw() +

  #eliminates background, gridlines, and chart border
  #theme(
   # plot.background = element_blank()
   #,panel.grid.major = element_blank()
   #,panel.grid.minor = element_blank()
  #)
p<-p+stat_smooth(method=lm, se=FALSE, color='black',alpha=0.1)
p<-p+theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20))
p <- p + theme(text = element_text(size = 20)) # this will change all text size 
#p<-p+theme(axis.title.x = element_text(size = rel(2)))+theme(axis.title.y = element_text(size = rel(2)))
p<-p + annotate('text', label=paste('R^2==',RSQ), parse=TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8) +annotate('text', label=paste('RMSE==',RMSE), parse=TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)


predicted.pq<-predict(rf.m,der1.ssn[,-1])
all.predicted<-cbind(all.predicted,predicted.pq)
fig.name=paste0(b,"/","Full_calibration_plots/",hd[q],".png")#This version does not save full calibration plots.
png(file=fig.name,height= 600,width=600)
p
dev.off()
}
#Combine the predicted values together
all.predicted.SSN<-cbind(as.vector(ssn),all.predicted)
colnames(all.predicted.SSN)<-c("SSN",hd)

colnames(msummary)<-c("Soil properties","OOB RMSEC","OOB Rsquared")
#Save full model summaries
write.table(msummary, file="Full models summary.csv",sep=",",row.names=FALSE)
#Save the linked file
write.table(all.predicted.SSN, file="All predictions.csv",sep=",", row.names=FALSE)
}
if(method=="PLS"){
	setwd(wd)
#download("https://www.dropbox.com/s/em8w2cdnnnkdrvo/rwanda.zip","rwanda.zip", mode="wb")
#unzip("rwanda.zip",overwrite=T)

#read raw IR data then preprocess by first derivatives using SG algorithm.
#mir<-read_csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/Calibration_Htsxt_MIR.csv")
mir<-infrared.data
ref<-reference.data
#Read reference data
#ref<-read_csv("~/Studies_data/data/Calibration_ref_data.csv")
#k<-which(mir$SSN%in%ref$ICRAF_ID)

#mir<-mir[k,]
#Step 2: Preprocess with SG algorithm
mir1<-as.matrix(mir[,-1])#Exclude metadata variables
wave<-as.numeric(substr(colnames(mir1),2,19))
colnames(mir1)<-wave
#First derivative
de1<-trans(mir1,tr="derivative",order=1,gap=23)
der1<-rev(as.data.frame(de1$trans))
colnames(der1)<-paste0("m",wave)

#Save derivative spectra
der1.ssn<-as.data.frame(cbind(as.vector(mir[,1]),der1)); colnames(der1.ssn)<-c("SSN",colnames(der1))
write.table(der1.ssn,file="First derivative.csv",sep=",",row.names=FALSE)
#merge with first derivative preprocessed spectra
ref.mir<-merge(ref,der1.ssn,by.x="SSN",by.y="SSN")
rc<-colnames(ref)
#which columns contains reference data?
ref<-ref.mir[,rc]
#Extract spectral predictors
mirp<-colnames(der1.ssn)[-1]
spectra<-der1.ssn[,mirp]
#Create two new subfolders within the current working using:
b<-getwd()
if(!file.exists("Models")){dir.create("Models")}
if(!file.exists("calibration_plots")){dir.create("calibration_plots")}

#Fit calibration models for the training set and use the testing set to validate the models.
#which are the training samples?
set.seed(67523)
#testing<-sample(1:nrow(ref.mir),0.33*nrow(ref.mir))#Hold-out a third of the calibration set for validation
#testing<-which(ref.mir$set=="val") #Use the val set rows
testing<-which(der1.ssn$SSN%in%hout$SSN)#with hout


#Use Kennard_Stone
sel <- kenStone(spectra,k=round(0.33*nrow(spectra)),pc=.99)
#View selected samples
plot(sel$pc[,1:2],xlab='PC1',ylab='PC2')
points(sel$pc[sel$model,1:2],pch=19,col=2) # points selected for calibration

#Loop for calibration of all soil properties in the reference set starts here
msummary<-NULL
hd<-colnames(ref)[-1]
for (q in 1:length(hd)){
refq<-which(colnames(ref)%in%hd[q])
ref.q<-ref[,refq]
pms.a <-NULL
pred.all <-NULL
cal<-cbind(as.vector(ref.q),spectra)[-testing,]
val<-cbind(as.vector(ref.q),spectra)[testing,]
colnames(cal)<-c(colnames(ref)[refq],colnames(spectra))
colnames(val)<-colnames(cal)
cal<-na.omit(cal)
val<-na.omit(val)
trainX <-cal[, -1]
trainY <-log(cal[,1])
tr0<-trainControl(method="LOOCV")
tr1<-trainControl(method = "repeatedcv", number = 10, repeats = 10, selectionFunction = "tolerance")
rf.m <- train(trainY~., method="pls", data=cal[,-1],trControl=tr1,tuneLength=20,metric="RMSE",preProc = c("center", "scale"))
	#Get final model to compute coefficient for variation explained
predi<-exp(predict(rf.m,rf.m$trainingData))
y<-cal[,1]
training.parameters<-round(postResample(predi,y),2)#computes RMSE and R-squared values for the calibration set
RSQ<-training.parameters[2]
RMSE<-training.parameters[1]
#Predict qth soil property of the holdoutset using  the MIR data and compare with the actual measurement
predi.test<-exp(predict(rf.m,val[,-1]))
y.test<-val[,1]
#Get PCs used
PCs<-rf.m$finalModel$ncomp
testing.parameters<-round(postResample(predi.test,y.test),2)#computes RMSE and R-squared values for the validation set
RSP<-testing.parameters[2]
RMSEP<-testing.parameters[1]
model.summary<-c(hd[q],PCs,training.parameters,testing.parameters)
msummary<-rbind(msummary,model.summary)
saveRDS(rf.m,file=paste0(b,"/","models/",hd[q],".rds"))
#Training
pm<-as.data.frame(cbind(y,predi))
colnames(pm)<-c("measured","predicted")
p<-ggplot(pm, aes(x=measured,y=predicted))+
geom_point(col="black",size=3,alpha=0.2)+
ggtitle(paste0("Calibration for ",hd[q]))+
	#xlim(as.numeric(as.vector(substr(variable,2,19))))+
	xlab("Measured")+
	ylab("Predicted")
	#theme with white background
  #theme_bw() +

  #eliminates background, gridlines, and chart border
  #theme(
   # plot.background = element_blank()
   #,panel.grid.major = element_blank()
   #,panel.grid.minor = element_blank()
  #)
p<-p+stat_smooth(method=lm, se=FALSE, color='black',alpha=0.1)
p<-p+theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20))
p <- p + theme(text = element_text(size = 20)) # this will change all text size 
#p<-p+theme(axis.title.x = element_text(size = rel(2)))+theme(axis.title.y = element_text(size = rel(2)))
p<-p + annotate('text', label=paste('R^2==',RSQ), parse=TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8) +annotate('text', label=paste('RMSE==',RMSE), parse=TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)

pmp<-as.data.frame(cbind(y.test,predi.test))#Validation data
colnames(pmp)<-c("measured.test","predicted.test")
p2<-ggplot(pmp, aes(x=measured.test,y=predicted.test))+
geom_point(col="brown",size=3,alpha=0.2)+
ggtitle(paste0("Validation for ",hd[q]))+
	#xlim(as.numeric(as.vector(substr(variable,2,19))))+
	xlab("Measured")+
	ylab("Predicted")
	#theme with white background
  #theme_bw() +

  #eliminates background, gridlines, and chart border
  #theme(
   # plot.background = element_blank()
   #,panel.grid.major = element_blank()
   #,panel.grid.minor = element_blank()
  #)
p2<-p2+stat_smooth(method=lm, se=FALSE, color='brown',alpha=0.1)
p2<-p2+theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20))
p2 <- p2 + theme(text = element_text(size = 20)) # this will change all text size 
#p<-p+theme(axis.title.x = element_text(size = rel(2)))+theme(axis.title.y = element_text(size = rel(2)))
p2<-p2 + annotate('text', label=paste('R^2==',RSP), parse=TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8) +annotate('text', label=paste('RMSE==',RMSEP), parse=TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)

png(file=paste0(b,"/Calibration_plots/",hd[q],".png"),height=400,width=800)
grid.arrange(p,p2,nrow=1)
dev.off()
}
colnames(msummary)<-c("Soil properties","PCs","LOOCV RMSEC","LOOCV Rsquared", "Holdout RMSEP","Holdout Rsquared")
write.table(msummary,file="Model_Summary.csv",sep=",",row.names=FALSE)

### All Samples #####
b<-getwd()
if(!file.exists("Full_Models")){dir.create("Full_Models")}
if(!file.exists("Full_calibration_plots")){dir.create("Full_calibration_plots")}
msummary<-NULL
hd<-colnames(ref[,-1])#Exclude SSN 
all.predicted<-NULL
for (q in 1:length(hd)) {
refq<-which(colnames(ref)%in%hd[q])
ref.q<-ref[,refq]
pms.a <-NULL
pred.all <-NULL
cal<-cbind(as.vector(ref.q),spectra)
colnames(cal)<-c(colnames(ref)[refq],colnames(spectra))
p<-which(is.na(der1.ssn[,1])==TRUE)
ifelse(length(p)>0,ssn<-der1.ssn[-p,1],ssn<-der1.ssn[,1])
ifelse(length(p)>0,der1.ssn<-der1.ssn[-p,],der1.ssn<-der1.ssn)
#Select training and testing sets
cal<-na.omit(cal)
trainX <-cal[, -1]
trainY <-log(cal[,1])
rf.m <- train(trainY~., method="pls", data=cal[,-1],trControl=tr1,tuneLength=20,metric="RMSE",preProc = c("center", "scale"))

	#Get final model to compute coefficient for variation explained
predi<-exp(predict(rf.m,rf.m$trainingData))
y<-cal[,1]
#Get PCs used
PCs<-rf.m$finalModel$ncomp

training.parameters<-c(hd[q],PCs,round(postResample(predi,y),3))
msummary<-rbind(msummary,training.parameters)
saveRDS(rf.m,file=paste0(b,"/","Full_Models/",hd[q],".rds"))
#Training
pm<-as.data.frame(cbind(y,predi))
colnames(pm)<-c("measured","predicted")
p<-ggplot(pm, aes(x=measured,y=predicted))+
geom_point(col="black",size=3,alpha=0.2)+
ggtitle(paste0("Calibration for ",hd[q]))+
	#xlim(as.numeric(as.vector(substr(variable,2,19))))+
	xlab("Measured")+
	ylab("Predicted")
	#theme with white background
  #theme_bw() +

  #eliminates background, gridlines, and chart border
  #theme(
   # plot.background = element_blank()
   #,panel.grid.major = element_blank()
   #,panel.grid.minor = element_blank()
  #)
p<-p+stat_smooth(method=lm, se=FALSE, color='black',alpha=0.1)
p<-p+theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=20))
p <- p + theme(text = element_text(size = 20)) # this will change all text size 
#p<-p+theme(axis.title.x = element_text(size = rel(2)))+theme(axis.title.y = element_text(size = rel(2)))
p<-p + annotate('text', label=paste('R^2==',RSQ), parse=TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8) +annotate('text', label=paste('RMSE==',RMSE), parse=TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)


predicted.pq<-predict(rf.m,der1.ssn[,-1])
all.predicted<-cbind(all.predicted,predicted.pq)
fig.name=paste0(b,"/","Full_calibration_plots/",hd[q],".png")#This version does not save full calibration plots.
png(file=fig.name,height= 600,width=600)
p
dev.off()
}
#Combine the predicted values together
all.predicted.SSN<-cbind(as.vector(ssn),all.predicted)
colnames(all.predicted.SSN)<-c("SSN",hd)

#Link the combined table to the sample field details.
 #if(file.exists("sample code.csv")){
#field<-read_csv("sample code.csv")
#field.predicted<-merge(field,all.predicted.SSN)}

 #if(!file.exists("sample code.csv")){
#field.predicted<-all.predicted.SSN}
colnames(msummary)<-c("Soil properties","PCs","repeatedCV RMSEC","repeatedCV Rsquared")
#Save full model summaries
write.table(msummary, file="Full models summary.csv",sep=",",row.names=FALSE)
#Save the linked file
write.table(field.predicted, file="All predictions.csv",sep=",", row.names=FALSE)
}
}