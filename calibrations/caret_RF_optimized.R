###  Script for doing calibration models using IR data
# Andrew Sila , April 2016

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

registerDoParallel()
getDoParWorkers()

### To install any required pacakge use:
#install.packages(c("package1, package2,...packagen"), dependencies=TRUE)

#Start by setting working directory
a<-"~/Studies_data/Mercy_calibrations"#Mac OSX
#a<-"D:/oaf"#Windows
setwd(a)
#download("https://www.dropbox.com/s/em8w2cdnnnkdrvo/rwanda.zip","rwanda.zip", mode="wb")
#unzip("rwanda.zip",overwrite=T)

#read raw IR data then preprocess by first derivatives using SG algorithm.
mir<-read_csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/Calibration_Htsxt_MIR.csv")
#Read reference data
ref<-read_csv("~/Studies_data/data/Calibration_ref_data.csv")
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
ref.mir<-merge(ref,der1.ssn,by.x="ICRAF_ID",by.y="SSN")
rc<-colnames(ref)
#which columns contains reference data?
ref<-ref.mir[,rc]
#Extract spectral predictors
mirp<-colnames(der1.ssn)[-1]
spectra<-ref.mir[,mirp]
#Create two new subfolders within the current working using:
b<-getwd()
if(!file.exists("Models")){dir.create("Models")}
if(!file.exists("calibration_plots")){dir.create("calibration_plots")}

#Fit calibration models for the training set and use the testing set to validate the models.
#which are the training samples?
set.seed(67523)
testing<-sample(1:nrow(ref.mir),0.33*nrow(ref.mir))#Hold-out a third of the calibration set for validation
testing<-which(ref.mir$set=="val") #Use the val set rows

#Use Kennard_Stone
sel <- kenStone(spectra,k=round(0.33*nrow(spectra)),pc=.99)
#View selected samples
plot(sel$pc[,1:2],xlab='PC1',ylab='PC2')
points(sel$pc[sel$model,1:2],pch=19,col=2) # points selected for calibration
testing<-sel$model

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
rf.m<-randomForest(trainY~.,data=cal[,-1],corr.bias=TRUE)
pred2 <- predict(rf.m ,cal)
plot(pred2,cal[,1])

pred3 <- predict(rf.m ,val)
plot(pred3,val[,1])

#rf.m <- train(trainY~., method="rf", data=cal[,-1],
			#metric="RMSE",
			#maximise=FALSE,
			#ntrees=3,# more tuning values requires to be tested
			#samplesize=0.5*nrow(cal),
			#corr.bias=TRUE,
			#trControl=trainControl(method="oob"),
			#tuneGrid=mgrid)
#Save the model in the currect working directory.
#trees<-as.numeric(rf.m$bestTune)[1]
#l<- which(rf.m$results$mtry==trees)

#Get final model to compute coefficient for variation explained
#predi<-predict(rf.m,rf.m$trainingData)
#y<-rf.m$finalModel$y
predi<-predict(rf.m ,cal)
y<-trainY
training.parameters<-round(postResample(predi,y),3)#computes RMSE and R-squared values for the calibration set
#Predict qth soil property of the holdoutset using  the MIR data and compare with the actual measurement
#predi.test<-predict(rf.m,val[,-1])
#y.test<-val[,1]
predi.test<-predict(rf.m ,val)
y.test<-val[,1]
testing.parameters<-round(postResample(predi.test,y.test),3)#computes RMSE and R-squared values for the validation set
model.summary<-c(hd[q],training.parameters,testing.parameters)
msummary<-rbind(msummary,model.summary)
saveRDS(rf.m,file=paste0(b,"/","models/",hd[q],".rds"))
#Training
labels1<-paste0("R^2 == ", format(training.parameters[2], digits=3))
labels2<-paste0("RMSEC = ",format(training.parameters[1], digits=3))
labels<-as.character(as.expression(c(labels1,labels2)))
#Get arbitrary positioning of plot lables
ab<-0.85*max(c(predi,y))
cb<-0.8*mean(range(y))
bb<-0.85*mean(range(y))
rfp<-as.data.frame(cbind(y,predi))

p1<- ggplot(data=rfp, aes(x = y, y =predi)) + 
	geom_smooth(method = "lm", se=FALSE, color = "grey", formula = y~x, lwd=0.3,alpha=0.5) +
	geom_point(color = "forestgreen",alpha=0.7,size=3)+
	ggtitle(paste("OOB: ",hd[q]," ; n =",nrow(rfp)))+
	xlab("Measured values")+
	ylab("Predicted values")+
	theme_bw()+
	geom_text(aes(label=as.character(labels1)),parse=TRUE,y=bb,x=ab,size=4,family="Times",data=data.frame())+
	geom_text(aes(label=as.character(labels2)),x=ab,y=cb,size=4,family="Times",data=data.frame())+

	theme(
	plot.background = element_blank()
	,panel.grid.major = element_blank()
	#,panel.grid.minor = element_blank()
	)
		p1<-p1+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
	#Vaidation set
labels12<-paste0("R^2 == ", format(testing.parameters[2], digits=3))
labels22<-paste0("RMSEP = ",format(testing.parameters[1], digits=3))
labels<-as.character(as.expression(c(labels12,labels22)))
ab<-0.85*max(c(predi.test,y.test))
cb<-0.8*mean(range(y.test))
bb<-0.85*mean(range(y.test))
rfp<-as.data.frame(cbind(y.test,predi.test))

p2<- ggplot(data=rfp, aes(x = y.test, y =predi.test)) + 
	geom_smooth(method = "lm", se=FALSE, color = "grey", formula = y~x, lwd=0.3,alpha=0.5) +
	geom_point(color = "brown",alpha=0.7,size=3)+
	ggtitle(paste("Validation: ",hd[q],"; n =",nrow(rfp)))+
	xlab("Measured values")+
	ylab("Predicted values")+
	theme_bw()+
	geom_text(aes(label=as.character(labels12)),parse=TRUE,y=bb,x=ab,family="Times",data=data.frame())+
	geom_text(aes(label=as.character(labels22)),x=ab,y=cb,family="Times",data=data.frame())+

	theme(
	plot.background = element_blank()
	,panel.grid.major = element_blank()
	#,panel.grid.minor = element_blank()
	)
	p2<-p2+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

source('~/local.github/sd4_codes/multiplot.R', chdir = TRUE)
png(file=paste0(b,"/Calibration_plots/",hd[q],".png"),height=400,width=800)
multiplot(p1,p2,cols=2)
dev.off()
}
colnames(msummary)<-c("Soil properties","OOB RMSEC","OOB Rsquared", "Holdout RMSEP","Holdout Rsquared")
write.table(msummary,file="Model_Summary.csv",sep=",",row.names=FALSE)

### All Samples #####
b<-getwd()
if(!file.exists("Full_Models")){dir.create("Full_Models")}
if(!file.exists("Full_calibration_plots")){dir.create("Full_calibration_plots")}
msummary<-NULL
hd<-colnames(ref[,-c(1:2)])#Exclude SSN 
all.predicted<-NULL
for (q in 1:length(hd)) {
refq<-which(colnames(ref)%in%hd[q])
ref.q<-ref[,refq]
pms.a <-NULL
pred.all <-NULL
cal<-cbind(as.vector(ref.q),spectra)
colnames(cal)<-c(colnames(ref)[refq],colnames(spectra))
#Select training and testing sets
cal<-na.omit(cal)
trainX <-cal[, -1]
trainY <-cal[,1]
#mgrid<-expand.grid(.mtry=seq(50,600,50))
rf.m<-randomForest(trainY~.,data=cal[,-1],corr.bias=TRUE)

#rf.m <- train(trainY~., method="rf", data=cal[,-1],
			#metric="RMSE",
			#maximise=FALSE,
			#ntrees=3,
			#corr.bias=TRUE,
			#trControl=trainControl(method="oob"),
			#tuneGrid=mgrid)
#Save the model in the currect working directory.
#trees<-as.numeric(rf.m$bestTune)[1]
#l<- which(rf.m$results$mtry==trees)
#Get final model to compute coefficient for variation explained
predi<-predict(rf.m,cal)
y<-trainY
training.parameters<-c(hd[q],round(postResample(predi,y),3))
msummary<-rbind(msummary,training.parameters)
saveRDS(rf.m,file=paste0(b,"/","Full_Models/",hd[q],".rds"))
#Training
labels1<-paste0("R^2 == ", format(training.parameters[2], digits=3))
labels2<-paste0("RMSEC = ",format(training.parameters[1], digits=3))
labels<-as.character(as.expression(c(labels1,labels2)))
ab<-0.75*max(c(predi,y))
bb<-0.4*mean(range(y))
cb<-0.3*mean(range(y))
rfp<-as.data.frame(cbind(y,predi))
p.all<- ggplot(data=rfp, aes(x = y, y =predi)) + 
	geom_smooth(method = "lm", se=FALSE, color = "grey", formula = y~x, lwd=0.3,alpha=0.5) +
	geom_point(color = "forestgreen",alpha=0.7,size=3)+
	ggtitle(paste("OOB: ",hd[q]))+
	xlab("Measured values")+
	ylab("Predicted values")+
	theme_bw()+
	geom_text(aes(label=as.character(labels1)),parse=TRUE,y=bb,x=ab,size=4,family="Times",data=data.frame())+
	geom_text(aes(label=as.character(labels2)),x=ab,y=cb,size=4,family="Times",data=data.frame())+

	theme(
	plot.background = element_blank()
	,panel.grid.major = element_blank()
	)
		p.all<-p.all+theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))     
#We need to predict the full data set using the full calibration data
predicted.pq<-predict(rf.m,cal)
all.predicted<-cbind(all.predicted,predicted.pq)
fig.name=paste0(b,"/","Full_calibration_plots/",hd[q],".png")#This version does not save full calibration plots.
png(file=fig.name,height= 600,width=600)
p.all
dev.off()
}
#Combine the predicted values together
all.predicted.SSN<-cbind(as.vector(der1.ssn[,1]),all.predicted)
colnames(all.predicted.SSN)<-c("SSN",hd)

#Link the combined table to the sample field details.
  if(file.exists("sample code.csv")){
field<-read_csv("sample code.csv")
field.predicted<-merge(field,all.predicted.SSN)}

  if(!file.exists("sample code.csv")){
field.predicted<-all.predicted.SSN}

colnames(msummary)<-c("Soil properties","OOB RMSEC","OOB Rsquared")
#Save full model summaries
write.table(msummary, file=paste(a,"/Full models summary.csv"),sep=",",row.names=FALSE)
#Save the linked file
write.table(field.predicted, file=paste0(a,"/All predictions.csv"),sep=",", row.names=FALSE)
