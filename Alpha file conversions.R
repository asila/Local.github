# Conversion of SPC files into a single flat table
#Andrew Sila
#A.sila@cgiar.org
#10.Oct.2011
###############################

setwd("/Volumes/NO NAME/test/chin");
tmpf<-list.files();
out<-"Alpha MIR.txt"
a<-which(tmpf=="LogFile.txt");
b<-which(tmpf==out);

if(length(a)!=0){tmpf<-tmpf[-a]};
if(length(b)!=0){tmpf<-tmpf[-b]};

waveb<-read.table("~/GRP4/R-scripts/SPC to filemaker new/Alpha filemaker column names.csv",sep=",",dec=".",header=F)[,1];

# Get the number of characters of each sample name:

test<-c();
for(i in 1:length(tmpf)){
	test[i]<-nchar(tmpf[i],type="chars")
	}
	
# Get the Batch_Labid from each sample name:

test.1<-c(rep(NA,length(tmpf)));
for(i in 1:length(tmpf)){
	test.1[i]<-substr(tmpf[i],1,(test[i]-6))
	}
	
# Create the additional information matrix:

add<-matrix(nrow=length(test.1),ncol=3,dimnames=list(test.1,c("Scan_date","Resolution","Zero_filling")));

library(hexView);

## Loop for reading spectra, making compatible and getting the additional information: 

m.c<-matrix(nrow=length(test.1),ncol=length(waveb),dimnames=list(test.1,waveb))
for(i in 1:length(test.1)){
	
# Read the .spc-file into R:

# Get the number of wavebands, the first and last waveband and the absorption values:
n<-readRaw(tmpf[i],offset=4,nbytes=4,human="int",size=4,endian="little")[[5]];
Ffirst<-readRaw(tmpf[i],offset=8,nbytes=8,human="real",size=8,endian="little")[[5]];
Flast<-readRaw(tmpf[i],offset=16,nbytes=8,human="real",size=8,endian="little")[[5]];
spec<-readRaw(tmpf[i],offset=544,nbytes=(4*n),human="real",size=4,endian="little")[[5]];
#if(Ffirst>Flast){c<-Ffirst;Ffirst<-Flast;Flast<-c;spec<-rev(spec);
#};
# Check if the wavebands have the same first waveband and number of wavebands compared to filemaker:
if(round(Ffirst,6)==waveb[1] & n==length(waveb))
# If yes, read the spectrum:
{	m.c[i,]<-spec;
}
# If the wavebands are different
else
{	# Read and calculate the differeing wavebands:
	
	bla<-c(rep(NA,n-2));
	for(k in 1:(n-2)){
		bla[k]<-k*(Flast-Ffirst)/(n-1);
		}
	w<-Ffirst+bla;
	wa<-round(c(Ffirst,w,Flast),6);
	# Make the spectrum compatible:
	m.c[i,]<-round(spline(wa,spec,method="natural",xout=waveb)[[2]],6);
	# Replace interpolations outside the measured range 		by missing values:
	for(l in 1:length(wa)){
	if((waveb[l]>wa[1])==F)break
	if(waveb[l]>wa[1])
	{
	m.c[i,l]<-NA	
	}}
		for(m in (length(waveb)):(length(waveb)-		length(wa))){
	if((waveb[m]<wa[length(wa)])==F)break;
	if(waveb[m]<wa[length(wa)])
	{
	m.c[i,m]<-NA	
	}}
}

test<-as.character(readRaw(tmpf[i],offset=32,nbytes=4,human="int",size=1.5,machine="binary",signed=F,endian="little"));
test.a<-paste(substr(test,35,42),substr(test,26,33),substr(test,17,21),sep="")
year<-substr(test.a,1,12);
a<-nchar(year)
b<-c(rep(NA,a));
k<-0;
for(j in a:1){
	b[j]<-as.numeric(substr(year,j,j))*2^k;
	k<-k+1
	}
year <-sum(b)
month<-substr(test.a,13,16);
a<-nchar(month)
b<-c(rep(NA,a));
k<-0;
for(j in a:1){
	b[j]<-as.numeric(substr(month,j,j))*2^k;
	k<-k+1
	}
month<-sum(b)
day<-substr(test.a,17,21);
a<-nchar(day)
b<-c(rep(NA,a));
k<-0;
for(j in a:1){
	b[j]<-as.numeric(substr(day,j,j))*2^k;
	k<-k+1
	}
day <-sum(b)
add[i,"Scan_date"]<-paste(month,"/",day,"/",year,sep="")

test<-readRaw(tmpf[i],offset=5830,nbytes= 360,human="char",endian="little");
test<-blockString(test);
test
test<-strsplit(test,"..",fixed=T);
add[i,"Resolution"]<-as.integer(substr(test[[1]][6],nchar(test[[1]])[6],nchar(test[[1]])[6]));
add[i,"Zero_filling"]<-as.integer(substr(test[[1]][18],nchar(test[[1]])[18],nchar(test[[1]])[18]));
	}

#Round off the decimal places in the column header to 1!
colnames(m.c)<-round(as.numeric(colnames(m.c)),1)

#Get laser wavenumber
laser.all<-c()
Laser_Wavelength<-c()
for (i in 1:length(tmpf)) {
test<-readRaw(tmpf[i],offset=11004,nbytes= 83,human="char",endian="little");

test<-blockString(test);test
test<-strsplit(test,"..",fixed=T);
lwn<-as.vector(substr(test[[1]],1,3))
y<-which(lwn=="LWN")
Laser_Wavelength[i]<-as.numeric(substr(test[[1]][y],6,nchar(test[[1]][y])))
laser.all<-c(laser.all,Laser_Wavelength[i])
	}

# Plot spectra:

plot(m.c[i,3:(ncol(m.c)-2)]~waveb[3:(ncol(m.c)-2)],type="l",xlim=c(max(waveb[3:(ncol(m.c)-2)]),min(waveb[3:(ncol(m.c)-2)])),ylim=c(min(m.c[,3:(ncol(m.c)-2)]),max(m.c[,3:(ncol(m.c)-2)])),col="purple");
for(i in 2:nrow(m.c)){lines(m.c[i,3:(ncol(m.c)-2)]~waveb[3:(ncol(m.c)-2)],col="purple")};

# Ask for user command y/n:

test<-menu(c("Spectra ok - please continue","Spectra not ok - please stop"),graphics=T);
graphics.off();
if(test==2)stop("You stopped importing those spectra");

#Get limits for plotting  the laser wavelength
icraf<-11604.1
daw<-Laser_Wavelength[1]

mn<-min(icraf,daw)
mx<-max(icraf,daw)
#Scaling for the laserwavelegth plot
a<-mn-5
b<-mx+5

#Check the Laserwavenumber
plot(laser.all,pch=19,cex=0.5,col="red",xlab="Samples scanned",ylab="Laser-Wavelength",main="Laser-wavelength trend",ylim=c(a,b), sub="Continous blue line=ICRAF Alpha")
abline(h=11604.1,col="blue",lty=1,lwd=2)


# Check the trend of the wavenumbers and decide:
test<-menu(c("Spectral Laser-Wavelength ok - please continue","Spectral Laser-Wavelength off - please tune up the instrument"),graphics=T)
graphics.off();
if(test==2)stop("Decide which spectra were completely off");


colnames(m.c)<-paste("a",colnames(m.c),sep="");
add.laser<-cbind(add,laser.all)
row.names(add.laser)<-tmpf
row.names(m.c)<-tmpf
colnames(add.laser)<-c(colnames(add),"Laser_Wavelength")
test<-merge(add.laser,round(m.c,6),by="row.names");

names(test)[1]<-c("SSN");

#Round-off decimal points in header
a<-which(substr(colnames(test),1,1)=="a")[1]

test$SSN<-substr(test$SSN,1,10) #Remove the  file extension .0spc
write.table(test,file=out,sep=",",dec=".",row.names=F);

