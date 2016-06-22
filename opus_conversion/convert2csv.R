#######################
##Converting opus files
#######################

## Read me
##This scripts works with spectral files named using a prefix 3-letter prefix followed
##by a six digit number

## 1. To use the script:
## 2.source the read.opus function
## 3.Change path to point to the folder with OPUS files
setwd("")#where read.opus.R file is saved
source('read.opus.R', chdir = TRUE)#This is the conversion program which needs to be sourced from where it is stored
opus.file.path<-''#where are opus files kept on your computer? Give full pathname
output.path<-''#give path name where converted file will be kept
file<-"raw spectra.csv"#Name to be given to converted OPUS into a single csv file

lst<- as.list(list.files(path=opus.file.path, pattern=".[0-9]$", full.names=TRUE))[-3170]
spectra<-c()
for ( i in 1:length(lst)){
spec <- opus(lst[[i]], speclib="ICRAF",plot.spectra=TRUE,print.progress=TRUE)
spectra<-rbind(spectra,spec)
}
#View part of spectra
spectra[,1:8]
setwd(output.path)
write.table(spectra,file=file,sep=",",row.names=FALSE)