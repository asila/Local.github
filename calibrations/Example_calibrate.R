source('~/local.github/calibrations/caret_RF_optimized_source.R', chdir = TRUE)
#1. Setwd
#2. Read MIR data
#3. Read reference data

wd<-"~/Studies_data/Mercy_calibrations/calibration_fn"
mir<-read_csv("~/Training/India/Support_Sinha/MIR_calval_data.csv")[,-c(2:7)]
colnames(mir)<-c("SSN",colnames(mir)[-1])
ref<-read_csv("~/Training/India/Support_Sinha/reference_calval_data.csv")[,-2]
calibrate(wd,mir,ref)