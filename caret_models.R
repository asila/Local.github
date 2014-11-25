#Load all the packages required
library(soil.spec)
library(pls)
library(chemometrics)
library(e1071)
library(caret)
library(spls)
#Source of models
#http://topepo.github.io/caret/modelList.html
#Start by specifying the folder to keep results. For me I chose "/Regional training/Nigeria/" and declare it as:

#Note for windows OS ensure the complete address path is defined inclusing the computer driver letter

results.folder<-"~/caret"

#Read into R both IR and chem data

raw<-read.table(file.choose(),sep=",",header=T)


nraw<-names(raw)
w2<-substr(colnames(raw[,200:250]),1,1)[20]

g<-which(nraw==paste(w2,2379.8,sep=""))
h<-which(nraw==paste(w2,2350.8,sep=""))

ifelse(length(g)>0,"The CO2 region 2379.8 to 2350.8 cm-1 will be removed","No CO2 bands found in this IR data")

ifelse(length(g)>0,raw<-raw[,-c(g:h)],raw<-raw)

cnr<-colnames(raw)
vari<-select.list(cnr,graphics=TRUE,multiple=TRUE,title="Columns independent variables...")
vari<-which(cnr==vari)
chem<-raw[,vari]

#spectra
cnr<-colnames(raw)
vari.ir<-select.list(cnr,graphics=TRUE,multiple=TRUE,title="Explanatory variables(IR)...")
spec<-raw[,which(cnr==vari.ir[1]):which(cnr==vari.ir[length(vari.ir)])]

#################################################################
#Create  folders for storing PLS plots and predicted results     #
	
dir.create(paste(results.folder,"PLS plots",sep="/"),showWarnings=FALSE)

dir.create(paste(results.folder,"PLS residual plots",sep="/"),showWarnings=FALSE)
	
dir.create(paste(results.folder,"Predictions",sep="/"),showWarnings=FALSE)

dir.create(paste(results.folder,"Calibration models",sep="/"),showWarnings=FALSE)


pls.d<-paste(results.folder,"PLS plots",sep="/")
resid.d<-paste(results.folder,"PLS residual plots",sep="/")
mvr.tn.d<-paste(results.folder,"Predictions",sep="/")
calib.d<-paste(results.folder,"Calibration models",sep="/")

#Set the spectra table as a matrix
mim<-as.matrix(spec)
dim(mim)
#Set the colnames of mim matrix as numeric
colnames(mim)<-as.numeric(substr(colnames(mim),2,15))

#msc preprocessing
raw.msc<-msc(mim)
dim(raw.msc)

#Combine msc with the ssn
mscssn<-cbind(as.vector(raw[,1]),raw.msc)
#View the first four columns and rows
mscssn[1:4,1:4]

#give label to the ssn column
colnames(mscssn)<-c("SSN",colnames(raw.msc))

mscssn<-as.data.frame(mscssn)
setwd(results.folder)
write.table(mscssn,"Mulplicative Scatter corrected spectra.csv",sep=",",row.names=F)


#Preprocess the raw spectra by first derivative; use library soil.spec
de<-trans(mim,tr="derivative",order=1,gap=21)
der<-as.matrix(de$trans)

dc<-ncol(der)

#Combine derivative with the ssn
derssn<-cbind(as.vector(raw[,1]),der)
#View the first four columns and rows
derssn[1:4,1:4]

#give label to the ssn column
colnames(derssn)<-c("SSN",colnames(der))

der<-as.data.frame(derssn)
setwd(results.folder)
write.table(derssn,"First derivative.csv",sep=",",row.names=F)

#Select either the first derivative or MSC file to use for calibration models
files<-c("First derivative.csv","Mulplicative Scatter corrected spectra.csv")
spe<-menu(files,graphics=TRUE,title="Select IR data")

#Read the file selected
der<-read.csv(paste(results.folder,files[spe],sep="/"),sep=",",header=T)

#First explore chem data and decide on the best transformation method

wtc<-colnames(chem)

#################################################
#Split the plot region into 4
par(mfrow=c(2,2))
plot.var<-menu(wtc,graphics=TRUE,title="Select variable to check distribution")
plotv<-as.vector(subset(chem,select=wtc[plot.var]))
hist(na.omit(plotv[-1,1]),main=wtc[plot.var])
plot(density(sqrt(na.omit(plotv[-1,1]))),main=wtc[plot.var])
################################################
cb<-cbind(chem,der)
cb<-cb[,-c(2:6)]
cfit<-train(Ca~.,data=cb,method="pls",tuneLength=10)
summary(cfit)
pred.m<-extractPrediction(list(cfit),testX=cb[,-1],testY=cb[,1])
pred.mt<-subset(pred.m,pred.m$dataType=="Training")
plotObsVsPred(pred.mt)

cb<-cbind(chem,der)
cb<-cb[,-c(1:2,4:6)]
cfit<-train(pH~.,data=cb,method="pls",tuneLength=10)
summary(cfit)
pred.m<-extractPrediction(list(cfit),testX=cb[,-1],testY=cb[,1])
pred.mt<-subset(pred.m,pred.m$dataType=="Training")
plotObsVsPred(pred.mt)

cb<-cbind(chem,der)
cb<-cb[,-c(1,3:6)]
cfit<-train(P~.,data=cb,method="pls",tuneLength=10)
summary(cfit)
pred.m<-extractPrediction(list(cfit),testX=cb[,-1],testY=cb[,1])
pred.mt<-subset(pred.m,pred.m$dataType=="Training")
plotObsVsPred(pred.mt)

 #cfit<-train(pH~.,data=cb,method="knn",tuneLength=10)
summary(cfit)
pred.m<-extractPrediction(list(cfit),testX=cb[,-1],testY=cb[,1])
plotObsVsPred(pred.m)

#################
cfit<-train(pH~.,data=cb,method="svmLinear")
summary(cfit)
pred.m<-extractPrediction(list(cfit),testX=cb[,-1],testY=cb[,1])
plotObsVsPred(pred.m)

cb<-calibset[,-c(1:4,6:20)]
#################
cfit<-train(K~.,data=cb,method="svmLinear",)
summary(cfit)
pred.m<-extractPrediction(list(cfit),testX=cb[,-1],testY=cb[,1])
plotObsVsPred(pred.m)

#################
cfit<-train(K~.,data=cb,method="earth")#knn; qrf is very poor
summary(cfit)
pred.m<-extractPrediction(list(cfit),testX=cb[,-1],testY=cb[,1])
plotObsVsPred(pred.m)