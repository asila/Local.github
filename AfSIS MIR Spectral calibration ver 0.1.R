#Script for creating PLS models using MIR data;
#Required files: MIR data and Reference data;
#Andrew Sila
#A.sila@cgiar.org
#21.March.2013
###############################

#Load all the packages required
library(soil.spec)
library(pls)
library(chemometrics)
library(e1071)

#Start by specifying the folder to keep results. For me I chose "~/GRP4/AfSIS-reporting/Versioned" and declare it as:

#Note for windows OS ensure the complete address path is defined inclusing the computer driver letter e.g. "D:/GRP4/AfSIS-reporting/Versioned"

results.folder<-"~/GRP4/Studies_data/Moi University/Results"

##################################################################
#Create  folders for storing PLS plots and predicted results     #
	
dir.create(paste(results.folder,"PLS plots",sep="/"),showWarnings=FALSE)
	
dir.create(paste(results.folder,"Predictions",sep="/"),showWarnings=FALSE)

dir.create(paste(results.folder,"Calibration models",sep="/"),showWarnings=FALSE)


pls.d<-paste(results.folder,"PLS plots",sep="/")
pred.d<-paste(results.folder,"Predictions",sep="/")
calib.d<-paste(results.folder,"Calibration models",sep="/")


#Read spectral data into raw;
raw<-read.csv("~/GRP4/Studies_data/Moi University/Raw MIR.csv.gz")

#Read reference data into chem;
chem<-read.csv("~/GRP4/Studies_data/Moi University/wet chemistry.csv.gz")

#Exclude the samples with SSNS icr030376, icr030377 and icr030532 which has high carbon values > 10%

re1<-which(chem$SSN==c("icr030376"))
re2<-which(chem$SSN==c("icr030377"))
re3<-which(chem$SSN==c("icr030532"))

#chem<-chem[-c(re1,re2,re3),]
nraw<-colnames(raw)

w2<-substr(colnames(raw[,200:250]),1,1)[20]
#Pick where spectra starts
a<-which(substr(nraw,1,1)==w2)[1]

#Ensure that columns for the spectra have been rounded of to one decimal point

numc<-paste(w2,round(as.numeric(substr(colnames((raw)[,-c(1:(a-1))]),2,15)),1),sep="")

#Set colnames to the rounded off
colnames(raw)<-c(colnames(raw)[1:(a-1)],numc)

nraww<-colnames(raw)
g<-which(nraww==paste(w2,2379.8,sep=""))
h<-which(nraww==paste(w2,2350.8,sep=""))

ifelse(length(g)>0,"The CO2 region 2379.8 to 2350.8 cm-1 will be removed","No CO2 bands found in this IR data")

ifelse(length(g)>0,raw<-raw[,-c(g:h)],raw<-raw)



#chem<-chem.50100
dim(chem)

#Remove the first column from raw spectra table which contains SSN/SSN and also the first and last columns of the spectra

rh<-colnames(raw)
rhw<-substr(rh,1,1)
l<-which(rhw==w2)[2]
z<-ncol(raw)-1

raw.cutoff<-raw[,l:z]

#Set the spectra table as a matrix
mim<-na.omit(as.matrix(raw.cutoff))
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
der<-read.table(paste(results.folder,files[spe],sep="/"),sep=",",header=T)
###----Merge first derivative with the chemical data

#First explore chem data and decide on the best transformation method

wtc<-colnames(chem)

#################################################
#Split the plot region into 4
par(mfrow=c(2,2))
plot.var<-menu(wtc,graphics=TRUE,title="Select variable to check distribution")
plotv<-as.vector(subset(chem,select=wtc[plot.var]))
hist(na.omit(plotv[-1,1]),main=wtc[plot.var])
plot(density(log(na.omit(plotv[-1,1]))),main=wtc[plot.var])
################################################
chem[1:4,1:4]
calibset<-merge(chem,der,by="SSN")
dim(calibset)

#If the two tables have the common field under different labels, we exclusively declare by.x and by.y
#calibset<-merge(chem,der,by.x="SSN",by.y="Batch_Labid")

#To know the cutoff regions get the columns of the calibset
colnames(calibset)

c<-ncol(chem)
wavenumbers<-as.numeric(substr(colnames(calibset[,-c(1:c)]),2,14))
y<-select.list(as.character(wavenumbers),graphics=TRUE,title="Select region where MIR begins and ends",multiple=TRUE)

ynum<-as.numeric(y)
ynumg<-max(ynum)
ynumn<-min(ynum)


a<-which(wavenumbers==ynumg)+c
b<-which(wavenumbers==ynumn)+c

#cutspec should lie in the regions 4000 to 600
cutspec<-as.matrix(calibset[,a:b])

y<-ncol(cutspec)

#Start PLS, by obtaining the variable names for cutchem matrix
#Add the library for calibration
#Check the cutoff IR regions
cutspec[1:4,1:4];cutspec[1:4,(y-4):y]

#Obtain from the calibset, the part with chem data
calibset[1:5,1:c]

#Where does the reference chem start?
refh<-colnames(calibset[,1:c])

z<-select.list(as.character(refh),graphics=TRUE,title="Select the first column where reference data starts",multiple=FALSE)
refh.z<-which(refh==z)

cutchem<-as.matrix(calibset[,refh.z:c])

colnames(cutspec)<-as.numeric(substr(colnames(cutspec),2,12))


#Create the prediction set from the first derivative table; ensuring the same regions in the cutspec are in the prediction set
der[1:4,1:4]

ps<-as.matrix(der[,(a-(c-1)):(b-(c-1))])

par(mfrow=c(1,1))
colnames(ps)<-colnames(cutspec)
dim(cutspec)

#############
#Natural log
############
summ.all2<-c()

for(q in 1:ncol(cutchem)){
setwd(pls.d)

hd<-colnames(cutchem)

#Include a vector for the units in the reference data and maintain the order
units<-c("%","%","%","%")
png(paste(hd[q],"calibration.png"),height=480,width=480,units="px")


cutspe<-as.data.frame(cutspec)
#Exclude constants
ifelse(min(na.omit(cutchem[,q]))==max(na.omit(cutchem[,q])),q<-q+1,q<-q)
#Determine the max pcs to be tested btn 1:20
cutchem[,q]<-ifelse(cutchem[,q]>0,cutchem[,q],NA)

pc<-c(1:20)
jj<-length(na.omit(cutchem[,q]))
	te<-c()
 for (u in 1:length(pc)) {
 	te[u]<-jj/pc[u]}
 
te.sel<-which(te>16)
p<-te.sel[length(te.sel)]

ifelse(jj<32,p<-10,p<-p)

cpe<-na.omit(cbind(as.vector(cutchem[,q]),cutspec))
colnames(cpe)<-c(hd[q],colnames(cutspec))
cpe[1:4,1:4]

ref<-as.numeric(as.matrix(cpe[,1]))
cutspe<-as.data.frame(cpe[,-1])
#Perform double cross-validation if calibration dataset<501 spectra
ifelse(jj<3500,res.pls<-mvr_dcv(log(ref)~.,data=cutspe,plot.opt=FALSE,ncomp=p,method="svdpc",repl=2,selstrat="diffnext",segments0=2,segment0.type="interleaved",na.action=na.omit,segments=2,segment.type="interleaved"),"Do nothing")

#Ensure we do not fit a model with less than 5 PCs.
af<-ifelse(res.pls$afinal<5,res.pls$afinal<-5,res.pls$afinal)


#Run the best model
mvr.tn<-plsr(na.omit(log(as.numeric(cutchem[,q]))~cutspec,ncomp=af,validation ="CV"))

nc<-af

rmse<-RMSEP(mvr.tn)$val[nc]

se<-round(rmse,2)

pred<-mvr.tn

#a<-min(exp(pred$fitted),na.omit(as.numeric(cutchem[,q])))

a<-round(min(na.omit(as.numeric(cutchem[,q])),na.omit(exp(pred$fitted.values[,,nc]))),1)

am<-round(max(na.omit(as.numeric(cutchem[,q])),na.omit(exp(pred$fitted.values[,,nc]))),1)

#am<-round(max(exp(pred$fitted),na.omit(as.numeric(cutchem[,q])),1))

b<-am+(0.05*am)
#b<-am

plot(na.omit(as.numeric(cutchem[,q]))~na.omit(exp(pred$fitted.values[,,nc])),pch=19,col="grey",ylab="",xlab="",main="",ylim=c(a,b),xlim=c(a,b))
bl<-lm(na.omit(as.numeric(cutchem[,q]))~na.omit(exp(pred$fitted.values[,,nc])))
abline(a=0,b=1,col="grey")
r<-round(summary(bl)$r.squared,2)

require(stats)
#r<-round(cor(exp(pred$fitted.values[,,nc]),na.omit(cutchem[,q]),method="pearson")^2,2)
mtext(paste(hd[q]," (using ",nc," PCs); ","n=",jj,sep=""),side=3,line=0.7,cex=1)
mtext(paste("Predicted ","(",units[q],")",sep=""),side=1,line=2.4,cex=0.8)
mtext(paste("Measured ","(",units[q],")",sep=""),side=2,line=2.4,cex=0.8)
legend("topleft",c(paste("r.squared",r,sep="="),paste("rmse=",se,sep="")),bty="n")

dev.off()

pds<-c()
summ2<-as.vector(c(paste(hd[q]," (",jj,")",sep=""),r,nc,se))
summ.all2<-rbind(summ.all2,summ2)
#To predict write the general syntax as: predict("model","new spectra"," optimal no. of PCs  selected in the model")
colnames(ps)<-colnames(cutspec)
pds<-predict(mvr.tn,ps,nc)
ssn<-as.matrix(der[,1])
pdss<-cbind(ssn,pds)
colnames(pdss)<-c("SSN",hd[q])

#Save loadings plot
#Create wd

setwd(pls.d)
lod<-as.matrix(mvr.tn$loadings)

wavl<-as.numeric(colnames(cutspec))

png(filename=paste(hd[q],"Loadings for PC ",nc,".png",sep=""),height=480,width=480,units="px")

ro<-round((nc/2+0.4),0)
ror<-round(sqrt(ro),0)+1
par(mfrow=c(ror,ror))

for (z in 1:nc){
plot(wavl,lod[,z],type="l",col=z,xlim=c(max(wavl),min(wavl)),ylab=paste("PC ",z," loading",sep=""),xlab=expression("Wavenumbers cm"^-1),main=hd[q])}
dev.off()
setwd(pred.d)
write.table(pdss,file=paste(hd[q],".csv"),sep=",",row.names=F)

#Save model into the calibration model folder for future predictions
setwd(calib.d)
save(mvr.tn,file=paste(hd[q],".rda",sep=""))

}



#Write model summary into a file
setwd(pls.d)
colnames(summ.all2)<-c("Property (n)","r-squared","PCs","rmse")
write.table(summ.all2,file="Model summary.csv",sep=",",row.names=F)


######################################################################
#Get the prediction and output in one file

setwd(pred.d)
pd<-list.files(pattern=".csv");
tt<-as.vector(pd)
ttt<-as.vector(strsplit(tt, " .csv"))
#tt<-substr(tt,11,40)
pdm<-read.table(pd[1],sep=",",header=T)
		pdm<-as.matrix(pdm[,1])
for(k in 1:length(pd)){
	pda<-read.table(pd[k],sep=",",header=T)
		pdm<-cbind(pdm,exp(pda[,2]))}
		colnames(pdm)<-c("SSN",ttt)
write.table(pdm,file="All edited predictions.csv",sep=",",row.names=F)

dim(pdm)
