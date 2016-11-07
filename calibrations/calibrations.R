library(readr)
cal<-read_csv("~/GRP4/Studies_data/Bernard Waruru/JNIRS_data/Thesis Cleaned calibration wetchem sets_converted_units.csv")
colnames(cal)
val<-read_csv("~/GRP4/Studies_data/Bernard Waruru/JNIRS_data/Thesis Cleaned valibration wetchem sets_converted_units.csv")
calp<-cal[,c("air-dry moisture%","pH 1:5 in water","total C %","total N %","Organic carbon %","extractable Fe mgkg-1","CEC-soil cmolkg-1","total clay%","water disp clay%","water disp sand %")]
Sama<-cal[,"macro mechanical%"]
cor(na.omit(cbind(Sama,calp[,4])))# this remains the same as previous so lets work with the calibration.s


######################################################
#With validation_holdout
######################################################

#0. Start by sourcing calibration file: caret_RF_optimized_source.R
#1. Set working directory
#2. Read MIR data
#3. Read reference data
require(readr)
require(dplyr)
#source('~/GRP4/OFRA/scripts/calibrations/caret_RF_optimized_source_hout.R', chdir = TRUE)
source('/Users/asila/GRP4/Studies_data/Bernard Waruru/JNIRS_scripts/RF_PLS_optimal.R', chdir = TRUE)
wd<-"~/GRP4/Studies_data/Bernard Waruru/JNIRS_MIR_PLS"
mir<-read_csv("~/GRP4/Studies_data/Bernard Waruru/JNIRS_data/Raw_MIR_spectra.csv")[,-c(2:1814)]
ref<-as.data.frame(cal[,c("SSN","air-dry moisture%","pH 1:5 in water","total C %","total N %","Organic carbon %","extractable Fe mgkg-1","CEC-soil cmolkg-1","total clay%","water disp clay%","water disp sand %")])
colnames(ref)<-c("SSN","mc","pH","TC","TN","OC","Fe","CEC","tClay","wdc","wds")
hout<-val[,c("SSN","air-dry moisture%","pH 1:5 in water","total C %","total N %","Organic carbon %","extractable Fe mgkg-1","CEC-soil cmolkg-1","total clay%","water disp clay%","water disp sand %")]
colnames(hout)<-colnames(ref)
colnames(mir)<-c("SSN",colnames(mir[,-1]))

ref<-rbind(ref,hout)
calibrate(wd,mir,ref,hout,method="PLS")
#NIR
wdn<-"~/GRP4/Studies_data/Bernard Waruru/JNIRS_NIR_PLS"
nir<-read_csv("~/GRP4/Studies_data/Bernard Waruru/JNIRS_data/Raw_NIR_spectra.csv")[,c(1,107:1143)]
p<-which(t(as.vector(nir[,1]))%in%t(as.vector(mir[,1])))
calibrate(wdn,nir[p,],ref,hout,method="PLS")


