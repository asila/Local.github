#Plan
###################

#1. read in data with raw spectra
#2. Identify spectra for the endmembers
#3. For each endmember identified compute pairwise angles between all samples.
#4.	Assign spectra to the same group with the endmember giving smallest angle.
#5. In different plots show the cosine angles between all samples per endmember 

library(lsa)
library(matrixcalc)
library(sna)
library(proxy)
library(lattice)
library(ggplot2)
tophumus<-read.csv("~andrewsila/Thesis/Data/Wkenya plus humus raw spectra.csv")
afsis<-read.csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/Calibration_Htsxt_MIR.csv")
pureclay<-read.csv("~/HTS-xt_spectra/B. Pure Clay/Averaged Pure_minerals.csv")
#Merge tophumus with afsis
colnames(tophumus)<-c("SSN",paste0("m",round(as.numeric(substr(colnames(tophumus[,-1]),2,19)),1)))
sm<-which(colnames(tophumus)%in%colnames(afsis))
#ir.all<-rbind(tophumus[,sm],afsis)

colnames(pureclay)<-c("SSN",paste0("m",round(as.numeric(substr(colnames(pureclay[,-1]),2,19)),1)))
smc<-which(colnames(pureclay)%in%colnames(afsis))

ir.all<-rbind(pureclay[,smc],tophumus[,sm],afsis)

#AfSIS with tophumus
afs<-ir.all[-c(1:165),]#also use in the SFI computation

#wkiemp with tophumus
wki<-ir.all[13:172,]

#Set afs to be ir ot wki
ir<-wki

#AfSIS with tophumus
test<-wki

#compute pairwise cosines for each spectra against all the spectra in test
endm<-ir.all[1:10,]
all.in.endm<-NULL
for ( i in 1:nrow(test)){
	sim.all<-NULL
for ( k in 1:nrow(endm)){
sim<-round(cosine(t(rbind(endm[k,-1],test[i,-1])))[1,2],8)
sim.all<-c(sim.all,sim)}
#mir.cos<-round(cosine(t(mir)),3)
endmember.class<-which(sim.all==max(sim.all))
i.end.member<-cbind(rep(test[i,1],length(endmember.class)),as.character(endm[endmember.class,1]))

all.in.endm<-rbind(all.in.endm,i.end.member)}
colnames(all.in.endm)<-c("SSN","endmember.ssn")
all.in.endm[,1]
setwd("~andrewsila/Thesis/cosine")
png(file="WKIEMP spectra with endmembers.png",height=600,width=600)
par(mfrow=c(2,2))
#1:Determine how many samples have been assigned in each endmember.
ssn<-pureclay[,1]
cls<-unique(all.in.endm[,2])
grps<-NULL
pure<-which(all.in.endm[,2]==cls[1])
test.1<-test[pure,]
wavenumbers<-as.numeric(substr(colnames(test.1[,-1]),2,19))
plot(wavenumbers,test.1[1,-1],type="l",ylim=c(min(ir.all[,-1]),max(ir.all[,-1])),main=paste0(cls[1],";n=",nrow(test.1)),col="red",xlim=c(max(wavenumbers),min(wavenumbers)),xlab=expression("Wavenumbers cm"^-1),ylab="Absorbance")
for ( l in 2:nrow(test.1)){
lines(wavenumbers,test.1[l,-1],col="red")}
clcc<-which(endm[,1]%in%cls[1])
lines(wavenumbers,endm[clcc,-1],col="blue",lwd=2)
legend("topright",c(as.character(cls[1]),"soil spectra"),col=c("blue","red"),lwd=c(2,1),bty="n")

#2:Determine how many samples have been assigned in each endmember.
cls<-unique(all.in.endm[,2])
grps<-NULL
pure<-which(all.in.endm[,2]==cls[2])
test.1<-test[pure,]
wavenumbers<-as.numeric(substr(colnames(test.1[,-1]),2,19))
plot(wavenumbers,test.1[1,-1],type="l",ylim=c(min(ir.all[,-1]),max(ir.all[,-1])),main=paste0(cls[2],";n=",nrow(test.1)),col="blue",xlim=c(max(wavenumbers),min(wavenumbers)),xlab=expression("Wavenumbers cm"^-1),ylab="Absorbance")
for ( l in 2:nrow(test.1)){
lines(wavenumbers,test.1[l,-1],col="blue")}
clcc<-which(endm[,1]%in%cls[2])
lines(wavenumbers,endm[clcc,-1],col="red",lwd=2)
legend("topright",c(as.character(cls[2]),"soil spectra"),col=c("red","blue"),lwd=c(2,1),bty="n")

#3:Determine how many samples have been assigned in each endmember.
cls<-unique(all.in.endm[,2])
grps<-NULL
pure<-which(all.in.endm[,2]==cls[3])

test.1<-test[pure,]
wavenumbers<-as.numeric(substr(colnames(test.1[,-1]),2,19))
plot(wavenumbers,test.1[1,-1],type="l",ylim=c(min(ir.all[,-1]),max(ir.all[,-1])),main=paste0(cls[3],";n=",nrow(test.1)),col="grey",xlim=c(max(wavenumbers),min(wavenumbers)),xlab=expression("Wavenumbers cm"^-1),ylab="Absorbance")
for ( l in 2:nrow(test.1)){
lines(wavenumbers,test.1[l,-1],col="grey")}
clcc<-which(endm[,1]%in%cls[3])
lines(wavenumbers,endm[clcc,-1],col="red",lwd=2)
legend("topright",c(as.character(cls[3]),"soil spectra"),col=c("red","grey"),lwd=c(2,1),bty="n")

#all
#all
#1:Determine how many samples have been assigned in each endmember.
cls<-unique(all.in.endm[,2])
grps<-NULL
pure<-which(all.in.endm[,2]==cls[1])
test.1<-test[pure,]
wavenumbers<-as.numeric(substr(colnames(test.1[,-1]),2,19))
plot(wavenumbers,test.1[1,-1],type="l",ylim=c(min(ir.all[,-1]),max(ir.all[,-1])),main=paste("Overlay all; n=",nrow(ir)),xlim=c(max(wavenumbers),min(wavenumbers)),xlab=expression("Wavenumbers cm"^-1),ylab="Absorbance")
for ( l in 2:nrow(test.1)){
lines(wavenumbers,test.1[l,-1],col="red")}
clcc<-which(endm[,1]%in%cls[1])
lines(wavenumbers,endm[clcc,-1],col="red",lwd=2)

#2:Determine how many samples have been assigned in each endmember.
cls<-unique(all.in.endm[,2])
grps<-NULL
pure<-which(all.in.endm[,2]==cls[2])

test.1<-test[pure,]
wavenumbers<-as.numeric(substr(colnames(test.1[,-1]),2,19))
#plot(wavenumbers,test.1[1,-1],type="l",ylim=c(min(ir.all[,-1]),max(ir.all[,-1])),main=paste0(cls[2],nrow(test.1)))
for ( l in 2:nrow(test.1)){
lines(wavenumbers,test.1[l,-1],col="blue")}
clcc<-which(endm[,1]%in%cls[2])
lines(wavenumbers,endm[clcc,-1],col="blue",lwd=2)
#3:Determine how many samples have been assigned in each endmember.
cls<-unique(all.in.endm[,2])
grps<-NULL
pure<-which(all.in.endm[,2]==cls[3])

test.1<-test[pure,]
wavenumbers<-as.numeric(substr(colnames(test.1[,-1]),2,19))
#plot(wavenumbers,test.1[1,-1],type="l",ylim=c(min(ir.all[,-1]),max(ir.all[,-1])),main=paste0(cls[3],nrow(test.1)))
for ( l in 2:nrow(test.1)){
lines(wavenumbers,test.1[l,-1],col="grey")}
clcc<-which(endm[,1]%in%cls[3])
lines(wavenumbers,endm[clcc,-1],col="grey",lwd=2)
legend("topright",c("mont","illite","halloy"),lty=1,col=c("red","blue","grey"),bty="n")
dev.off()

#########################################################################################################
#Afsis
######################################
test<-afs
#compute pairwise cosines for each spectra against all the spectra in test
endm<-ir.all[1:10,]
all.in.endm<-NULL
for ( i in 1:nrow(test)){
	sim.all<-NULL
for ( k in 1:nrow(endm)){
sim<-round(cosine(t(rbind(endm[k,-1],test[i,-1])))[1,2],8)
sim.all<-c(sim.all,sim)}
#mir.cos<-round(cosine(t(mir)),3)
endmember.class<-which(sim.all==max(sim.all))
i.end.member<-cbind(rep(test[i,1],length(endmember.class)),as.character(endm[endmember.class,1]))

all.in.endm<-rbind(all.in.endm,i.end.member)}
colnames(all.in.endm)<-c("SSN","endmember.ssn")

#Extract details for the mineralogical classes and then store them.
ssni<-as.numeric(all.in.endm[,1])
all.in.endm->spec.grps
spec.grps[,1]<-as.vector(ir.all[,1][ssni])
spec.grps<-as.data.frame(spec.grps)
#Set working directory
setwd("~andrewsila/Thesis/cosine")
write.table(spec.grps,file="AfSIS spec classes.csv",sep=",",row.names=FALSE)
png(file="AfSIS spectra with endmembers.png",height=600,width=600)
par(mfrow=c(2,2))
#1:Determine how many samples have been assigned in each endmember.
ssn<-pureclay[,1]
cls<-unique(all.in.endm[,2])
grps<-NULL
pure<-which(all.in.endm[,2]==cls[1])
test.1<-test[pure,]
wavenumbers<-as.numeric(substr(colnames(test.1[,-1]),2,19))
plot(wavenumbers,test.1[1,-1],type="l",ylim=c(min(ir.all[,-1]),max(ir.all[,-1])),main=paste0(cls[1],";n=",nrow(test.1)),col="red",xlim=c(max(wavenumbers),min(wavenumbers)),xlab=expression("Wavenumbers cm"^-1),ylab="Absorbance")
for ( l in 2:nrow(test.1)){
lines(wavenumbers,test.1[l,-1],col="red")}
clcc<-which(endm[,1]%in%cls[1])
lines(wavenumbers,endm[clcc,-1],col="blue",lwd=2)
legend("topright",c(as.character(cls[1]),"soil spectra"),col=c("blue","red"),lwd=c(2,1),bty="n")

#2:Determine how many samples have been assigned in each endmember.
cls<-unique(all.in.endm[,2])
grps<-NULL
pure<-which(all.in.endm[,2]==cls[2])
test.1<-test[pure,]
wavenumbers<-as.numeric(substr(colnames(test.1[,-1]),2,19))
plot(wavenumbers,test.1[1,-1],type="l",ylim=c(min(ir.all[,-1]),max(ir.all[,-1])),main=paste0(cls[2],";n=",nrow(test.1)),col="blue",xlim=c(max(wavenumbers),min(wavenumbers)),xlab=expression("Wavenumbers cm"^-1),ylab="Absorbance")
for ( l in 2:nrow(test.1)){
lines(wavenumbers,test.1[l,-1],col="blue")}
clcc<-which(endm[,1]%in%cls[2])
lines(wavenumbers,endm[clcc,-1],col="red",lwd=2)
legend("topright",c(as.character(cls[2]),"soil spectra"),col=c("red","blue"),lwd=c(2,1),bty="n")

#3:Determine how many samples have been assigned in each endmember.
cls<-unique(all.in.endm[,2])
grps<-NULL
pure<-which(all.in.endm[,2]==cls[3])

test.1<-test[pure,]
wavenumbers<-as.numeric(substr(colnames(test.1[,-1]),2,19))
plot(wavenumbers,test.1[1,-1],type="l",ylim=c(min(ir.all[,-1]),max(ir.all[,-1])),main=paste0(cls[3],";n=",nrow(test.1)),col="grey",xlim=c(max(wavenumbers),min(wavenumbers)),xlab=expression("Wavenumbers cm"^-1),ylab="Absorbance")
for ( l in 2:nrow(test.1)){
lines(wavenumbers,test.1[l,-1],col="grey")}
clcc<-which(endm[,1]%in%cls[3])
lines(wavenumbers,endm[clcc,-1],col="red",lwd=2)
legend("topright",c(as.character(cls[3]),"soil spectra"),col=c("red","grey"),lwd=c(2,1),bty="n")

#4:Determine how many samples have been assigned in each endmember.
cls<-unique(all.in.endm[,2])
grps<-NULL
pure<-which(all.in.endm[,2]==cls[4])
test.1<-test[pure,]
wavenumbers<-as.numeric(substr(colnames(test.1[,-1]),2,19))
plot(wavenumbers,test.1[1,-1],type="l",ylim=c(min(ir.all[,-1]),max(ir.all[,-1])),main=paste0(cls[4],";n=",nrow(test.1)),col="grey",xlim=c(max(wavenumbers),min(wavenumbers)),xlab=expression("Wavenumbers cm"^-1),ylab="Absorbance")
for ( l in 2:nrow(test.1)){
lines(wavenumbers,test.1[l,-1],col="grey")}
clcc<-which(endm[,1]%in%cls[4])
lines(wavenumbers,endm[clcc,-1],col="red",lwd=2)
legend("topright",c(as.character(cls[4]),"soil spectra"),col=c("red","grey"),lwd=c(2,1),bty="n")

dev.off()

png("pure mineral spectra.png")
par(mfrow=c(3,4))
for (i in 1:nrow(endm)){
plot(wavenumbers,endm[i,-1],type="l",ylim=c(min(endm[,-1]),max(endm[,-1])),main=paste0(endm[i,1]),col=i,xlim=c(max(wavenumbers),min(wavenumbers)),xlab=expression("Wavenumbers cm"^-1),ylab="Absorbance")	
}
dev.off()