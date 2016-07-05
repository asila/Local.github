require("prospectr")
spectra<-read_csv("~/Training/India/Erick_Towett/raw spectra.csv")[,-c(1:7)]
#Use Kennard_Stone
png("~/Thesis/Plots/Kennard_demo.png",height=900,width=900)
par(mfrow=c(2,2))
sep=c(2,3,4,30)
for (i in 1:length(sep)){
sel <- kenStone(spectra,k=sep[i],pc=0.76)
#View selected samples
plot(sel$pc[,1:2],xlab='PC1',ylab='PC2',pch=19,col="blue",main=paste0(sep[i]," extreme points selected"))
points(sel$pc[sel$model,1:2],pch=19,col="red") # points selected for calibration
testing<-sel$model
}
dev.off()