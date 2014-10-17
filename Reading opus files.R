require(soil.spec)
require(chemometrics)

files<-"/Volumes/HTS_XT_AFSIS/from_ciesin/Hts-xt/ICRAF/RAW/MIR REFLECTANCE"#Rescans
lst.fo<-list.files(files)

spectra<-c()
meta<-c()
for ( k in 1:length(lst.fo)){
setwd(paste(files,lst.fo[k],sep="/"))

lstr<-as.list(list.files(pattern=".[0-9]$"))#With repeats
xx <- read.opus(lstr,plot.spectra=TRUE,speclib="ICRAF",print.progress=TRUE) 
spectra<-rbind(spectra,xx@data@ab)
meta<-rbind(meta,xx@data@samples)
}

#cbind spectra with meta
spectra<-cbind(meta,spectra[,-1])