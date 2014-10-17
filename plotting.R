mir<-read.csv(file.choose())

dim(mir)

mir[1:6,1:10]

wavenumber<-colnames(mir[,6:2547])
wavenumber
wavenumbers<-as.numeric(substr(wavenumber,2,7))

spectra<-as.matrix(mir[,6:2547])#or
class(spectra)

spectra[1:3,1:3]

tiff("Kari.tiff")
plot(wavenumbers,spectra[1,],type="l",ylab="Absorbance",col=2)

for ( i in 1:nrow(spectra)){
lines(wavenumbers,spectra[i,],col=i)}
dev.off()

