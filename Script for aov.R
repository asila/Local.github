#A simple method for doing ANOVA
cg<-read.csv(file.choose())
str(cg)
#Remove column 1
cg<-cg[,-1]
#Make ST, REP, ML and TAR as factors
cg$ST<-as.factor(cg$ST)
cg$REP<-as.factor(cg$REP)
cg$ML<-as.factor(cg$ML)
cg$TAR<-as.factor(cg$TAR)

cg[1:8,]
ch4<-aov(CH4~ST+ML+TAR,data=cg)
summary(ch4)

co2<-aov(CO2~ST*ML*TAR,data=cg)
summary(co2)

#Tonsee the means use
model.tables(co2)

#To test the means use:
TukeyHSD(co2,"ST")