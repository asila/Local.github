#Read HTS-xt
mir<-read.csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/Calibration_Htsxt_MIR.csv")
nir<-read.csv("~/Dropbox/AfSIS_reporting_data/Seperated_datasets/Calibration_MPA_NIR.csv")
dim(mir)

p<-which(nir[,1]%in%mir[,1])
nirp<-nir[p,-c(2:5)]


view.nir<-function(ir){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[1,-1],type="l",ylim=c(0,3),ylab="Absorbance",main="NIR spectra",pch="",xlim=c(8000,4000))
for ( i in 2:nrow(ir)){
	lines(wavenumbers,ir[i,-1],col="red")
}}
view.mir<-function(ir){
wavenumbers<-as.numeric(substr(colnames(ir[,-1]),2,19))
plot(wavenumbers,ir[1,-1],type="l",ylim=c(0,3),ylab="Absorbance",main="MIR spectra",pch="",xlim=c(4000,600))
for ( i in 2:nrow(ir)){
	lines(wavenumbers,ir[i,-1],col="red")
}}


nir1<-nir #remove non-spectra
waves<-as.numeric(substr(colnames(nir1[,-1]),2,19))
n<-which(waves>4000&waves<8000)
nir2<-nir[,c(1,n)]

waves<-as.numeric(substr(colnames(mir[,-1]),2,19))
m<-which(waves>600&waves<4000)
mir1<-mir[,c(1,m)]

nirmir<-merge(rev(nir2),mir1,by.x="SSN",by.y="SSN")
wavenumbers<-as.numeric(substr(colnames(nirmir[,-1]),2,19))

view.nirmir<-function(nirmir){
wavenumbers<-as.numeric(substr(colnames(nirmir[,-1]),2,19))
plot(wavenumbers,nirmir[1,-1],type="l",ylim=c(0,3),ylab="Absorbance",main="MIR spectra",pch="",xlim=c(8000,600))
for ( i in 2:nrow(nirmir)){
	lines(wavenumbers,nirmir[i,-1],col="grey")
}}
par(mfrow=c(1,3))
view.nir(nirp)
view.mir(mir)
view.nirmir(nirmir)
write.table(nirmir,file="~/Thesis/Data/AfSIS Fused_nir_mir.csv",sep=",",row.names=FALSE)


der[1:5,1:5]
dp<-prcomp(der[,-1])
plot(dp$x)
dpx<-as.data.frame(cbind(as.vector(der[,1]),dp$x))
colnames(dpx)<-c("SSN",colnames(dp$x))
codes<-read.csv(file.choose())
dpc<-merge(codes,dpx)
plot(dpc$PC1,dpc$PC2)
text(dpc$PC1,dpc$PC2,labels=dpc$Site)

write.table(dpc,file="afsis_PCA.csv",sep=",",row.names=FALSE)
dpc<-read.csv("afsis_PCA.csv")

plot(dpc$PC1,dpc$PC2,pch="")
sites<-unique(dpc$Site)
for ( i in 1:length(sites)){
	dps<-subset(dpc,dpc$Site==sites[i])
	text(dps$PC1,dps$PC2,col=i,labels=sites[i])
}
text(dpc$PC1,dpc$PC3,labels=dpc$Site)

summary(dp)