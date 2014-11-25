library(lsa)
library(matrixcalc)
library(sna)
library(proxy)
library(lattice)
library(ggplot2)
#source("http://www.phaget4.org/R/myImagePlot.R")
tophumus<-read.csv("~Andrewsila/Thesis/Data/Wkenya plus humus plus white sand raw spectra.csv")
afsis<-read.csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/Calibration_Htsxt_MIR.csv")
xrd<-read.csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/AfSIS_xrd.csv")
xrdm<-which(xrd$Montmorillonite==max(na.omit(xrd$Montmorillonite)))
xrdk<-which(xrd$Kaolinite==max(na.omit(xrd$Kaolinite)))

#Merge tophumus with afsis
colnames(tophumus)<-c("SSN",paste0("m",round(as.numeric(substr(colnames(tophumus[,-1]),2,19)),1)))
sm<-which(colnames(tophumus)%in%colnames(afsis))
ir.all<-rbind(tophumus[,sm],afsis)
#AfSIS with tophumus
afs<-ir.all[-c(1:160),]

#wkiemp with tophumus
wki<-ir.all[1:166,]

#Set afs to be ir ot wki
ir<-wki

l<-ncol(ir)
#ir<-ir[,-c(2:11,l)]
mir<-rev(ir[,-1])#Using lsa; ensure no missing values
mir.cos<-round(cosine(t(mir)),3)
#mir.cos<-round(dist(mir,method="manhattan"),3)

# get lower triangle using sna library
low<-upper.triangle(mir.cos)
low[1:16,1:16]
zero<-which(low==0)
loww<-low[-zero]
#Get minimum distance
min.c<-min(loww)

#Get max distance
max.c<-max(loww)

#Determine which pairs of samples has minumum dist
sep.spec<-which(low==min.c)
rss.css.all<-c()
for (k in 1:length(sep.spec)){
	css<-(sep.spec[k]/ nrow(low))#determine which column
	if(round(css,0)>sep.spec[k]%% nrow(low))
	{css<-round(css,0)+1}
	
	if(!round(css,0)>sep.spec[k]%% nrow(low))
	{css<-round(css,0)+1
	
	#find from the cssth column the corresponding row with simialrity distance value = min.c
	rss<-which(low[,css]==min.c)
	rss.css<-c(rss,css)
	rss.css.all<-list(rss.css.all,rss.css)}
	}
	css<-(sep.spec[k]/ nrow(low))#determine which column

	if(round(css,0)==nrow(low))
	{css<-round(css,0)
		rss<-sep.spec[k]%% nrow(low)
		rss.css.all<-c(rss,css)}


#Get unique pairs
#rss.u<-unique(rss.css.all[[-1]])
rss.u<-unique(rss.css.all)

#ir<-mir
#Plot the pairs
close.screen(all.screens=TRUE)
split.screen(c(1,2),erase=TRUE)
split.screen(c(2,2),erase=TRUE)

if(length(rss.u)>2){
	ros<-rss.u[-length(rss.u)]
	cols<-rss.u[length(rss.u)]
	view.spectra<-function(ir){
	wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
	plot(wavenumbers,ir[1,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra for most different pair(s)",col="blue")
		for ( i in 2:(nrow(ir)-1)){
		lines(wavenumbers,ir[i,-1],col="blue")
		}
		lines(wavenumbers,ir[nrow(ir),-1],col="red")
	}
screen(3)
irr<-which(ir[,1]=="")
view.spectra(ir[rss.u,])
}

if(!length(rss.u)>2){
	ros<-rss.u[1]
	cols<-rss.u[2]
	view.spectra<-function(ir){
	wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
	plot(wavenumbers,ir[1,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra for most different pair(s)",col="blue")
		
		lines(wavenumbers,ir[nrow(ir),-1],col="red")
	}
screen(3)

view.spectra(ir[rss.u,])
}

if(length(rss.u)<2){
	ros<-rss.u[1]
	cols<-rss.u[2]
	view.spectra<-function(ir){
	wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
	plot(wavenumbers,ir[1,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra for most different pair(s)",col="blue")
		
	}
screen(3)

view.spectra(ir[rss.u,])
}

legend("topleft",c("Pure humus","Very extreme spectra from pure humus;SL"),lty=1,col=c("red","blue"),bty="n")
setwd("~andrewsila/Thesis/cosine/wkenya")
png("wkiemp_pure humus plus extreme.png")
view.spectra(ir[rss.u,])
legend("topleft",c("Pure humus; SH","white sand; SL"),lty=1,col=c("red","blue"),bty="n")
dev.off()
#First the spectra picked under the column with min dist
if(!css==nrow(low)){
dist1<-c(low[1:(css-1),rss.u[k]],low[css,(css+1):nrow(low)])}
if(css==nrow(low)){
dist1<-c(low[1:(css-1),rss.u[2]])}

sf1<-as.data.frame(cbind(as.vector(ir[-css,1]),dist1))
if(!rss.u[1]==nrow(low)){
dist2<-c(low[1:(rss.u[1]-1),rss.u[1]],low[rss.u[1],(rss.u[1]+1):nrow(low)])}
if(rss.u[1]==nrow(low)){
dist2<-c(low[1:(rss.u[1]-1),rss.u[1]])}

sf2<-as.data.frame(cbind(as.vector(ir[-rss[1],1]),dist2))
sf12<-merge(sf1,sf2)

write.table(sf12,file="sfi.csv",sep=",",row.names=FALSE)
sf12<-read.csv("sfi.csv")
g1<-which(sf12$dist1>sf12$dist2)
g2<-which(!sf12$dist1>sf12$dist2)

#Get all neighbours for each of the pairs in most distant spectra
#par(mfrow=c(1,1))
#plot(sf12[,-1])

screen(4)
if(length(g1)<1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[css,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
}
if(length(g1)>1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[css,-1],pch="n",type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
for ( i in 1:length(g1)){
	lines(wavenumbers,ir[g1[i],-1],col="red")}
	lines(wavenumbers,ir[css,-1],col="black")

}

#g2
screen(5)
if(length(g2)<1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[rss.u[1],-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
}
if(length(g2)>1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[rss.u[1],-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="",lty=3)
for ( i in 1:length(g2)){
	lines(wavenumbers,ir[g2[i],-1],col="red")}
	lines(wavenumbers,ir[rss.u[1],-1],col="black")
}

screen(6)
if(length(g2)<1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[rss.u[1],-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
}
if(length(g2)>1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[rss.u[1],-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="",lty=3)
for ( i in 1:length(g2)){
	lines(wavenumbers,ir[g2[i],-1],col="red")}
	lines(wavenumbers,ir[rss.u[1],-1],col="black")
}

if(length(g1)<1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
#plot(wavenumbers,ir[css,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
}
if(length(g1)>1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
#plot(wavenumbers,ir[css,-1],pch="n",type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
for ( i in 1:length(g1)){
	lines(wavenumbers,ir[g1[i],-1],col="red")}
	lines(wavenumbers,ir[css,-1],col="black")

}

screen(2)
plot(sf12[,2],sf12[,3],pch="",xlab="Extreme s1 cosine distance", ylab="Extreme s2 cosine distance")
points(sf12[g1,2],sf12[g1,3],pch=16,col="red")
points(sf12[g2,2],sf12[g2,3],pch=16,col="blue")
legend("bottomleft",c("Angle cos for s1 and close points","Angle cosine for s2 and close points"),pch=19,col=c("red","blue"),bty="n")
ids<-identify(sf12[,2],sf12[,3],labels=ir[-css,1])

cod<-qplot(sf12[,2],sf12[,3],colour=sf12[,3],ylim=c(min(sf12[,3]),1),xlim=c(min(sf12[,2]),1),xlab="Cosine of angle between white sand (SH) and others",ylab="Cosine of angle between pure humus and others",main="Similarity measure of WKIEMP MIR spectra and pure humus")
cod + scale_colour_gradient(limits=c(min(sf12[,3]),1),low="red",high="blue",name="SFI")

png("SFI_wkiemp_pure humus plus ws.png")
cod + scale_colour_gradient(limits=c(min(sf12[,3]),1),low="red",high="blue",name="SFI")
dev.off()

test.samples<-ir[-rss.u,1]
sfi<-as.data.frame(cbind(as.vector(test.samples),grp1.dist,grp2.dist))
names(sfi)<-c("SSN","SFI1","SFI2")

lab<-read.csv("~/Thesis/Data/Spatial/labdata.csv")
lab.sfi<-merge(lab,sfi,by="SSN")
write.table(lab.sfi,file="~/Thesis/Data/labdata sfi.csv",sep=",",row.names=FALSE)


lab.sfi<-read.csv("~/Thesis/Data/labdata sfi.csv")
plot(lab.sfi[,c(32,31)])
######################################################################################################################
#AfSIS
######################################################################################################################
#Set afs to be ir ot wki

setwd("~andrewsila/Thesis/cosine/afsis")
#from endmember computation
afsis<-ir.all[-c(13:172),]
#ir<-afs
ir<-afsis
l<-ncol(ir)
#ir<-ir[,-c(2:11,l)]
mir<-rev(ir[,-1])#Using lsa; ensure no missing values
mir.cos<-round(cosine(t(mir)),3)
colnames(mir.cos)<-afsis[,1]
row.names(mir.cos)<-afsis[,1]
write.table(mir.cos,file="Cosine angle for Afsis samples.csv",sep=",",row.names=FALSE)
mir.cos<-read.csv("Cosine angle for Afsis samples.csv")
# get lower triangle using sna library
low<-upper.triangle(as.matrix(mir.cos))
low[1:16,1:16]
#calculate SFI1 and SFI2
s1<-which(colnames(mir.cos)=="icr119455")
s2<-which(colnames(mir.cos)=="Quartz")
sfi<-as.data.frame(cbind(as.vector(colnames(mir.cos)),mir.cos[,s1],mir.cos[,s2]))[-c(s1,s2),]
colnames(sfi)<-c("SSN","SFI","SFI2")
write.table(sfi,file="AfSiS humus quartz sfi.csv",sep=",",row.names=FALSE)
sfi<-read.csv("AfSiS humus quartz sfi.csv")
plot(sfi[,-1])
ide<-identify(sfi[,-1])

plot(sfi[-ide,-1])
zero<-which(low==0)
loww<-low[-zero]
#Get minimum distance
min.c<-min(loww)

#Get max distance
max.c<-max(loww)

#Determine which pairs of samples has minumum dist
sep.spec<-which(low==min.c)
rss.css.all<-c()
for (k in 1:length(sep.spec)){
	css<-(sep.spec[k]/ nrow(low))#determine which column
	if(round(css,0)>sep.spec[k]%% nrow(low))
	{css<-round(css,0)+1}
	#find from the cssth column the corresponding row with simialrity distance value = min.c
	rss<-which(low[,css]==min.c)
	rss.css<-c(rss,css)
	rss.css.all<-list(rss.css.all,rss.css)
	
		if(!round(css,0)>sep.spec[k]%% nrow(low))
	{css<-round(css,0)+1

	#find from the cssth column the corresponding row with simialrity distance value = min.c
	rss<-which(low[,css]==min.c)
	rss.css<-c(rss,css)
	rss.css.all<-list(rss.css.all,rss.css)}
	}
	
	#Get unique pairs
#rss.u<-unique(rss.css.all[[-1]])
rss.u<-as.vector(unique(rss.css.all)[[2]])

#ir<-mir
#Plot the pairs
close.screen(all.screens=TRUE)
split.screen(c(1,2),erase=TRUE)
split.screen(c(2,2),erase=TRUE)

if(length(rss.u)>2){
	ros<-rss.u[-length(rss.u)]
	cols<-rss.u[length(rss.u)]
	view.spectra<-function(ir){
	wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
	plot(wavenumbers,ir[1,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra for most different pair(s)",col="blue")
		for ( i in 2:(nrow(ir)-1)){
		lines(wavenumbers,ir[i,-1],col="blue")
		}
		lines(wavenumbers,ir[nrow(ir),-1],col="red")
	}
screen(3)
irr<-which(ir[,1]=="")
view.spectra(ir[rss.u,])
}

if(!length(rss.u)>2){
	ros<-rss.u[1]
	cols<-rss.u[2]
	view.spectra<-function(ir){
	wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
	plot(wavenumbers,ir[1,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra for most different pair(s)",col="blue")
		
		lines(wavenumbers,ir[nrow(ir),-1],col="red")
	}
screen(3)

view.spectra(ir[rss.u,])
}

if(length(rss.u)<2){
	ros<-rss.u[1]
	cols<-rss.u[2]
	view.spectra<-function(ir){
	wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
	plot(wavenumbers,ir[1,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra for most different pair(s)",col="blue")
		
	}
screen(3)

view.spectra(ir[rss.u,])
}

legend("topleft",c("Pure humus","white sand;SL"),lty=1,col=c("red","blue"),bty="n")
setwd("~andrewsila/Thesis/cosine/afsis")
png("afsis_pure humus plus extreme ws.png")
view.spectra(ir[rss.u,])
legend("topleft",c("AfSIS_Sanza Pombo sample","white sand;SL"),lty=1,col=c("red","blue"),bty="n")
dev.off()
#First the spectra picked under the column with min dist
if(!css==nrow(low)){
dist1<-as.data.frame(cbind(as.vector(ir[-1541,1]),c(mir.cos[1:1540,1541],mir.cos[1542:nrow(mir.cos),1541])))
dist2<-as.data.frame(cbind(as.vector(ir[-6,1]),c(mir.cos[1:5,6],mir.cos[7:nrow(mir.cos),6])))

sf12<-merge(dist1,dist2,by="V1")
colnames(sf12)<-c("SSN","sf1","sf2")
write.table(sf12,file="sfi.csv",sep=",",row.names=FALSE)
sf12<-read.csv("sfi.csv")
g1<-which(sf12$sf1>sf12$sf2)
g2<-which(!sf12$sf1>sf12$sf2)

#Get all neighbours for each of the pairs in most distant spectra
#par(mfrow=c(1,1))
#plot(sf12[,-1])

screen(4)
if(length(g1)<1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[css,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
}
if(length(g1)>1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[css,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
irs<-ir[-c(css,rss.u),]#Remove the extreme spectra
for ( i in 1:length(g1)){
	lines(wavenumbers,irs[g1[i],-1],col="red")}
	lines(wavenumbers,irs[css,-1],col="black")

}

#g2
screen(5)
if(length(g2)<1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[rss.u[1],-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
}
if(length(g2)>1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[rss.u[1],-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="",lty=3)
for ( i in 1:length(g2)){
	lines(wavenumbers,irs[g2[i],-1],col="red")}
	lines(wavenumbers,irs[rss.u[1],-1],col="black")
}

screen(6)
if(length(g2)<1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[rss.u[1],-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
}
if(length(g2)>1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[rss.u[1],-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="",lty=3)
for ( i in 1:length(g2)){
	lines(wavenumbers,ir[g2[i],-1],col="red")}
	lines(wavenumbers,ir[rss.u[1],-1],col="black")
}

if(length(g1)<1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
#plot(wavenumbers,ir[css,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
}
if(length(g1)>1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
#plot(wavenumbers,ir[css,-1],pch="n",type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
for ( i in 1:length(g1)){
	lines(wavenumbers,ir[g1[i],-1],col="red")}
	lines(wavenumbers,ir[css,-1],col="black")

}

screen(2)
plot(sf12[,2],sf12[,3],pch="",xlab="Extreme s1 cosine distance", ylab="Extreme s2 cosine distance")
points(sf12[g1,2],sf12[g1,3],pch=16,col="red")
points(sf12[g2,2],sf12[g2,3],pch=16,col="blue")
legend("bottomleft",c("Angle cos for s1 and close points","Angle cosine for s2 and close points"),pch=19,col=c("red","blue"),bty="n")
ids<-identify(sf12[,2],sf12[,3],labels=ir[-css,1])

cod<-qplot(sf12[,2],sf12[,3],colour=sf12[,2],ylim=c(min(sf12[,3]),0.8),xlim=c(min(sf12[,2]),1),xlab="Cosine of angle between Sanza_Pombo sample and others",ylab="Cosine of angle between WS and others",main="Similarity measure AfSIS MIR spectra and white sand")
png("AfSIS-ws SFI.png")
cod + scale_colour_gradient(limits=c(min(sf12[,2]),1),low="red",high="blue",name="SFI")
dev.off()
############################################################################################
#fix pure earth

########################### AfSIS samples closeness to high carbon black earth
#use black earth sample with highest carbon
p<-which(ir$SSN=="icr119455")
mp<-min(low[p,-c(1:(1+p))])
q<-which(low[p,]==mp)[[1]]

css<-p
rss.u<-c(q,p)

#ir<-mir
#Plot the pairs
close.screen(all.screens=TRUE)
split.screen(c(1,2),erase=TRUE)
split.screen(c(2,2),erase=TRUE)

if(length(rss.u)>2){
	ros<-rss.u[-length(rss.u)]
	cols<-rss.u[length(rss.u)]
	view.spectra<-function(ir){
	wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
	plot(wavenumbers,ir[1,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra for most different pair(s)",col="blue")
		for ( i in 2:(nrow(ir)-1)){
		lines(wavenumbers,ir[i,-1],col="blue")
		}
		lines(wavenumbers,ir[nrow(ir),-1],col="red")
	}
screen(3)
irr<-which(ir[,1]=="")
view.spectra(ir[rss.u,])
}

if(!length(rss.u)>2){
	ros<-rss.u[1]
	cols<-rss.u[2]
	view.spectra<-function(ir){
	wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
	plot(wavenumbers,ir[1,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra for most different pair(s)",col="blue")
		
		lines(wavenumbers,ir[nrow(ir),-1],col="red")
	}
screen(3)

view.spectra(ir[rss.u,])
}

if(length(rss.u)<2){
	ros<-rss.u[1]
	cols<-rss.u[2]
	view.spectra<-function(ir){
	wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
	plot(wavenumbers,ir[1,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra for most different pair(s)",col="blue")
		
	}
screen(3)

view.spectra(ir[rss.u,])
}

legend("topleft",c("Pure humus","white sand;SL"),lty=1,col=c("red","blue"),bty="n")
setwd("~andrewsila/Thesis/cosine/afsis")
png("afsis_pure humus plus ws.png")
view.spectra(ir[rss.u,])
legend("topleft",c("Pure humus","white sand;SL"),lty=1,col=c("red","blue"),bty="n")
dev.off()
#First the spectra picked under the column with min dist
if(!css==nrow(low)){
dist1<-as.data.frame(cbind(as.vector(ir[-4,1]),c(mir.cos[1:3,4],mir.cos[5:nrow(mir.cos),4])))
dist2<-as.data.frame(cbind(as.vector(ir[-6,1]),c(mir.cos[1:5,6],mir.cos[7:nrow(mir.cos),6])))

sf12<-merge(dist1,dist2,by="V1")
colnames(sf12)<-c("SSN","sfi1","sfi2")
write.table(sf12,file="sfi2.csv",sep=",",row.names=FALSE)
sf12<-read.csv("sfi2.csv")
g1<-which(sf12$sfi1>sf12$sfi2)
g2<-which(!sf12$sfi1>sf12$sfi2)

#Get all neighbours for each of the pairs in most distant spectra
#par(mfrow=c(1,1))
#plot(sf12[,-1])

screen(4)
if(length(g1)<1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[css,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
}
if(length(g1)>1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[css,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
irs<-ir[-c(css,rss.u),]#Remove the extreme spectra
for ( i in 1:length(g1)){
	lines(wavenumbers,irs[g1[i],-1],col="red")}
	lines(wavenumbers,irs[css,-1],col="black")

}

#g2
screen(5)
if(length(g2)<1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[rss.u[1],-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
}
if(length(g2)>1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[rss.u[1],-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="",lty=3)
for ( i in 1:length(g2)){
	lines(wavenumbers,irs[g2[i],-1],col="red")}
	lines(wavenumbers,irs[rss.u[1],-1],col="black")
}

screen(6)
if(length(g2)<1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[rss.u[1],-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
}
if(length(g2)>1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[rss.u[1],-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="",lty=3)
for ( i in 1:length(g2)){
	lines(wavenumbers,ir[g2[i],-1],col="red")}
	lines(wavenumbers,ir[rss.u[1],-1],col="black")
}

if(length(g1)<1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
#plot(wavenumbers,ir[css,-1],type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
}
if(length(g1)>1){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
#plot(wavenumbers,ir[css,-1],pch="n",type="l",ylim=c(min(ir[,-1]),max(ir[,-1])),xlim=c(max(wavenumbers),min(wavenumbers)),ylab="Absorbance",main="Spectra close to extreme s1",pch="")
for ( i in 1:length(g1)){
	lines(wavenumbers,ir[g1[i],-1],col="red")}
	lines(wavenumbers,ir[css,-1],col="black")

}

screen(2)
plot(sf12[,2],sf12[,3],pch="",xlab="Extreme s1 cosine distance", ylab="Extreme s2 cosine distance")
points(sf12[g1,2],sf12[g1,3],pch=16,col="red")
points(sf12[g2,2],sf12[g2,3],pch=16,col="blue")
legend("bottomleft",c("Angle cos for s1 and close points","Angle cosine for s2 and close points"),pch=19,col=c("red","blue"),bty="n")
ids<-identify(sf12[,2],sf12[,3],labels=ir[-css,1])

cod<-qplot(sf12[,2],sf12[,3],colour=sf12[,2],ylim=c(min(sf12[,3]),0.8),xlim=c(min(sf12[,2]),1),xlab="Cosine of angle between pure humus(SH) and others",ylab="Cosine of angle between WS and others",main="Similarity measure AfSIS MIR spectra and pure humus")
 cod + scale_colour_gradient(limits=c(min(sf12[,2]),1),low="red",high="blue",name="SFI")
 png("tophumus-ws SFI.png")
cod + scale_colour_gradient(limits=c(min(sf12[,2]),1),low="red",high="blue",name="SFI")
dev.off()

 afsis.lab<-read.csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/AfSIS reference data.csv")
lab.sfi<-merge(afsis.lab,sf12,by="SSN")
write.table(lab.sfi,file="~/Thesis/Data/AfSIS_labdata sfi.csv",sep=",",row.names=FALSE)

plot(lab.sfi$sfi1,lab.sfi$sfi2)
lab.sfi<-read.csv("~/Thesis/Data/AfSIS_labdata sfi.csv")
plot(lab.sfi[,c(58,10)])
png("Carbon vs SFI.png")
plot(lab.sfi[,c(58,28)],xlab="SFI")
dev.off()
png("m3.al vs SFI.png")
plot(lab.sfi[,c(58,13)],xlab="SFI")
dev.off()

names(lab.sfi)

library(FactoMineR)
lab.sfi<-read.csv("~/Thesis/Data/AfSIS_labdata sfi.csv")

labs<-lab.sfi[,-c(1:11)]
#Read XRD data
xrd<-read.csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/AfSIS_xrd.csv")
lbx<-merge(lab.sfi,xrd)

labs<-lbx[,-c(1:11)]
minerals<-c("Quartz","Kaolinite","Microcline","Albite","Montmorillonite")
min<-which(colnames(labs)%in%minerals)
sfi.pca<-PCA(labs[,c(1:17,32:34,47,min)],scale.unit=TRUE,ind.sup=labs[,1])
labs$ssfi1