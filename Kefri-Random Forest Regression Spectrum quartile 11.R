# RANDOM FOREST REGRESSION
# Keith D Shepherd 2012


## STEPS: (1) Set working directory. (2) Select variables. (3) Run and review the model. (4) Reduce no of trees and no of variables. (5) Rerun the model. (6) Name version and save. (7) Print results of interest to file. 

## Note: normally essential code sections are marked with ##; required user inputs are flagged with ####; and optional sections with #.

## Load required libraries
library(soil.spec)
library(randomForest)
library(ggplot2)
library(grid)
library(fields)
#library(rpanel)
library(Hmisc)
library(reshape)
library(gtools)
library(rattle)
library(latticist)
#library(DiagnosisMed) (tcl screen pops up which freezes makes R to stop running; run this library at the diagnostics step, towards the end of this script )
library(gregmisc)

#### Set working directory
setwd("~/GRP4/Studies_data/Kefri/RandomForests")
dp<-getwd()

##################################################################
#Create  folders for storing Random forest graphics and the calibration model     #
	
dir.create(paste(dp,"Diagnostics plots",sep="/"),showWarnings=FALSE)
	
dir.create(paste(dp,"Models",sep="/"),showWarnings=FALSE)


dia.d<-paste(dp,"Diagnostics plots",sep="/")
mod.d<-paste(dp,"Models",sep="/")

## Clean out old dataframes and variables
rm(list=ls());



#### Specify comma-delimted data file
data.inr <- read.table("~/Dropbox/AfSIS reporting data/23May2012_Hts-xt MIR calibration.csv.gz", sep=",",header=T)
data.in <- with(data.inr,{data.inr[!(Country=="Haiti"),]}) #

#Read data from IR table and reference and merge
ir<-read.csv(file.choose())
renw<-paste("m",round(as.numeric(substr(colnames(ir[,-1]),2,17)),1),sep="")
colnames(ir)<-c("SSN",renw)
ref<-read.csv(file.choose())
data.in<-merge(ref,ir,by="SSN")
#data.in <- with(data.inr,{data.inr[(Study=="AfSIS core"),]}) # Use instead of above for AfSIS core only
#write.table(data.in, file=paste("datacheck.csv",sep="") , sep=",",row.names=F)

#### Any additional variable calcs
data.in <- within(data.in,{
	#socind <- (Acidified.Carbon - (0.1 + 1.0 * wcvfrfc10)) / ((1.0 + 10 * wcvfrfc10) - (0.1 + 1.0 * wcvfrfc10))
	ExAl <- m3.Al/90
	Alsat <- ExAl/(ExBas+ExAl)
	Acsat <- ExAc/(ExAc+ExBas)
})

#### Specifiy and check spectral data and soil data; check for missing values in spectral data. Ensure to inlcude calculated soil vars. Version 1
# colnames(data.in) 
data.spec <- subset(data.in, select = c(m4001.6:m601.7)) # Omit bands with NA
data.soil <- subset(data.in, select = c(SSN:PSI,ExAc:Acsat))

#Select the variables interactively- Version 2
dataspec <- as.data.frame(select.list(as.character(colnames(data.in)), graphics=TRUE,title="IR region",multiple=TRUE)) # select the first and last bands,omit intial bands with NA

datasoil <- select.list(as.character(colnames(data.in)), graphics=TRUE,title="Exclude IR",multiple=TRUE) # Omit bands with NA

#Get spectral data and soil property data matrices
data.spec <- data.in[,which(colnames(data.in)==dataspec[1,1]):which(colnames(data.in)==dataspec[2,1])]# Omit bands with NA
data.soil <- data.in[,-c(which(colnames(data.in)==dataspec[1,1]):which(colnames(data.in)==dataspec[2,1]))]

# Checks to soil data
soil.stats <- as.data.frame(t(stats(data.soil)))
soil.stats <- subset(soil.stats, select = c(N,min,max))
soil.stats
checkmiss <- as.data.frame(t(stats(data.spec)))
checkmiss[which(checkmiss[,9] > 0),]
datacheck <- data.soil[complete.cases(data.soil),]
# datacheck
soil.stats <- as.data.frame(t(stats(datacheck)))
soil.stats <- subset(soil.stats, select = c(N,min,max))
soil.stats # Complete cases only

#Preprocess the raw spectra by the selected method
detach(package:gregmisc)
detach(package:gplots)

colnames(data.spec) <- as.numeric(substr(colnames(data.spec),2,15)) # Convert wavebands to numerical
de<-trans(data.spec, tr="derivative" , order=1 , gap=21) # When using continuum removal the parameters for order and gap are ignored.
der<-as.matrix(de$trans) # Note trans relabels wavebands with X
data.trans<-data.frame(data.soil,der) #Combine derivative with the soil data

# Row numbers check
nrow(der)
nrow(data.soil)
nrow(data.trans)
#data.trans <- within(data.trans,{soci <- (Total.Carbon/(100-psa.c4sand))}) #calculate new variables here
#data.trans <- subset(data.trans, soci<0.5,) # subset to remove outliers
#plot(data.trans$psa.c4sand,data.trans$wcvfrairdry) #Plot variables of interest
#subset(data.in,Site == "Kisongo", select = c(a.brooks:ksat)) # Examine input data for a site


#### Assign target variable and predictor variables #############################
#yy <- "ExBas" # define y variable

# Specify spectral variables as range and non-spectral variables by name without quotes or column no. Default m601.7:m4001.6 for 1st Derivs and m603.6:m3997.8 for continuum removed.

# Variable list 
#Depth,Cultivated,#pH, m3.Al, m3.B, m3.Cu ,m3.Fe, m3.Mn, m3.P, m3.S, m3.Zn, PSI,#ExNa, ExCa, ExMg, ExK, ExBas, ESR, ESP, ECd, ExAc, CaMg,#llwcgpct,lshrinkpct,piwcgpct,plwcgpct,#a.brooks,b.brooks,alpha.brooks,wcvfri,sucjkgi,wcvfrsat,wcvfrfc10,wcvfrfc33,wcvfrwp1500,wcvfrairdry,awc1,awc2,ksat,#psa.asand,psa.asilt,psa.aclay,#psa.c1sand,psa.c1silt,psa.c1clay,psa.c2sand,psa.c2silt,psa.c2clay,psa.c3sand,psa.c3silt,psa.c3clay,psa.c4sand,psa.c4silt,psa.c4clay,#psa.w1sand,psa.w1silt,psa.w1clay,psa.w2sand,psa.w2silt,psa.w2clay,psa.w3sand,psa.w3silt,psa.w3clay,psa.w4sand,psa.w4silt,psa.w4clay,#Total.Nitrogen,Total.Carbon,Acidified.Nitrogen,Acidified.Carbon,#Na,Mg,Al,P,S,Cl,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,As,#Br,Rb,Sr,Y,Zr,Ba,La,Ce,Pr,Nd,Sm,Hf,Ta,W,Pb,Bi,Th,#Quartz,Albite,Microcline,Hornblende,Tremolite,Diopside,muscovite,Montmorillonite,Clinochlore,#Kaolinite,Halloysite,Vermiculite,Hematite,Goethite,Magnetite,Gibbsite,Ilmenite,Anatase,Calcite,Gypsum,Dolomite,#X601.7:X4001.6

### Select the X variables for the forest model.  EXCLUDE THE Y VARIABLE

#yyy<-c("pH","m3.Al","m3.B","m3.Cu","m3.Fe","m3.Mn","m3.P","m3.S","m3.Zn","PSI","ExNa","ExCa","ExMg","ExK","ExBas","ESR","ESP","ECd","ExAc","CaMg","llwcgpct","lshrinkpct","piwcgpct","plwcgpct","a.brooks","b.brooks","alpha.brooks","wcvfri","sucjkgi","wcvfrsat","wcvfrfc10","wcvfrfc33","wcvfrwp1500","wcvfrairdry","awc1","awc2","ksat","psa.asand","psa.asilt","psa.aclay","psa.c1sand","psa.c1silt","psa.c1clay","psa.c2sand","psa.c2silt","psa.c2clay","psa.c3sand","psa.c3silt","psa.c3clay","psa.c4sand","psa.c4silt","psa.c4clay","psa.w1sand","psa.w1silt","psa.w1clay","psa.w2sand","psa.w2silt","psa.w2clay","psa.w3sand","psa.w3silt","psa.w3clay","psa.w4sand","psa.w4silt","psa.w4clay","Total.Nitrogen","Total.Carbon","Acidified.Nitrogen","Acidified.Carbon","Na","Mg","Al","P","S","Cl","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","As","Br","Rb","Sr","Y","Zr","Ba","La","Ce","Pr","Nd","Sm","Hf","Ta","W","Pb","Bi","Th","Quartz","Albite","Microcline","Hornblende","Tremolite","Diopside","muscovite","Montmorillonite","Clinochlore","Kaolinite","Halloysite","Vermiculite","Hematite","Goethite","Magnetite","Gibbsite","Ilmenite","Anatase","Calcite","Gypsum","Dolomite")
yyy<-c("pH","m3.Al","m3.B","m3.Cu","m3.Fe","m3.Mn","m3.P","m3.S","m3.Zn","PSI","ExNa","ExCa","ExMg","ExK","ExBas","ESR","ESP","ECd","ExAc","CaMg")

#setwd("~/AfSIS-reporting/ Random forest models")

for ( j in 1:length(yyy)){
	
yy<-yyy[j]
#yy<-"Total.Carbon" ; if devloping one randomforest at a time, set the name of the property here.
	
	a<-as.data.frame(paste("X",colnames(der)[1],sep=""))[1,1]
	b<-as.data.frame(paste("X",colnames(der)[ncol(der)],sep=""))[1,1]
	
	aa<-which(colnames(data.trans)==a)
	bb<-which(colnames(data.trans)==b)
	
	datax <- data.trans[,aa:bb]
	
	################################################################################
	
	## Prepare desired data and retain only complete cases
	datat1 <- data.frame(data.trans[yy],datax) # Combine Y and X variables with Y in the first column
	datatr <- datat1[complete.cases(datat1),] # retain only complete cases
	nrow(datatr) # show number of rows of complete data
	
	## Prepare metadata file for target variable
	meta = data.in[complete.cases(data.in[,yy]),]
	
	#soil.starts <- select.list(as.character(colnames(meta)), graphics=TRUE,title="Soil data starts...",multiple=TRUE) # Omit bands with NA
	
	#Get spectral data and soil property data matrices
	#meta.ends <- which(colnames(meta)==soil.starts)-1
	
	#metadat <- meta[,1:meta.ends]
	#nrow(metadat)
	#nrow(datatr)
	
	## Run Random Forest model. Note: it is better not to use model formulation style on large data sets.
	set.seed(131)
	datatr.rf <- randomForest(data=datatr, x=datatr[,-1], y=datatr[,1], ntree=200, importance=TRUE, na.action=na.omit)
	datatr.rf
	
	varImpPlot(datatr.rf)
	
	# Plot OOB residuals
	length(datatr.rf$y)
	length(datatr.rf$predicted)
	nrow(metadat)
	resids <- data.frame(metadat, datatr.rf$y, datatr.rf$predicted)
	resids <- rename(resids, c(datatr.rf.y="actual")) 
	resids <- rename(resids, c(datatr.rf.predicted="predicted")) 
	resids <- within(resids, {residuals <- (predicted - actual)})
	plot(resids$actual, resids$residuals)
	
	#### Panel plot by site.
	ggplot(resids,aes(x=actual, y=residuals, colour = Depth)) + geom_point() + xlab("Actual") + ylab("Residuals") + facet_wrap(~ Site, ncol = 6)
	
	# Interactive plot
	latticist(resids)
	
	## Plot predictions using the model
	modres1 <- predict(datatr.rf, datatr)
	model.pred<-paste("Model predicted ",yy)
	oob<-paste("Out of Bag ",yy)
	
	mp1 <- ggplot(datatr, aes(modres1, datatr[,yy])) + geom_point(shape=1) + geom_smooth(method = "lm") + xlab("Predicted") + ylab("Actual") +
	  opts(title =model.pred,  axis.text.x = theme_text(size=11),
	  axis.title.x = theme_text(hjust=0.55, vjust = 0, size=12),
	  axis.text.y = theme_text(size=11),
	  axis.title.y = theme_text(hjust=0.6, vjust=0.35,angle=90,size=12),
	  plot.title=theme_text(hjust=0.55,vjust=1,size=13,face="bold"),
	  panel.border=theme_rect(size=0.75),
	  legend.position="none"
	   )
	# Plot predictions using out-of-bag
	mv1 <- ggplot(datatr, aes(datatr.rf$predicted, datatr.rf$y)) + geom_point(shape=1) + geom_smooth(method = "lm") + xlab("Predicted") + ylab("Actual") + 
	  opts(title =oob, axis.text.x = theme_text(size=11),
	  axis.title.x = theme_text(hjust=0.55, vjust = 0, size=12),
	  axis.text.y = theme_text(size=11),
	  axis.title.y = theme_text(hjust=0.6, vjust=0.35,angle=90,size=12),
	  plot.title=theme_text(hjust=0.55,vjust=1,size=13,face="bold"),
	  panel.border=theme_rect(size=0.75),
	  legend.position="none"
	   )
	# Print plots together
	grid.newpage()
	vport <- function(x, y) viewport(layout.pos.row=x,layout.pos.col=y);
	pushViewport(viewport(layout = grid.layout(nrow=1, ncol=2, widths=unit(3.5,"inches"), heights=unit(3.5,"inches"))))
	print(mp1, vp=vport(1,1));
	print(mv1, vp=vport(1,2));
	
	## Print statistics of actual vs predicted fits for training model
	mp <- lm(datatr[,yy] ~ modres1) # Linear regression stats
	mp.train <- list(summary(mp)$adj.r.squared, summary(mp)$sigma, summary(mp)$coef)
	mp.train
	#calculate the RMSE for actual vs predicted for training model
	(mean((datatr[,yy] - modres1)^2))^0.5
	
	## Print statistics of actual vs OOB prediction
	mpoob <- lm(datatr.rf$y ~ datatr.rf$predicted) # Linear regression stats
	mp.oob <- list(summary(mpoob)$adj.r.squared, summary(mpoob)$sigma, summary(mpoob)$coef)
	mp.oob
	#calculate the RMSE for actual vs predicted OOB
	(mean((datatr.rf$y - datatr.rf$predicted)^2))^0.5
	
	
	## Organize and plot"importance" of variables: higher means more important: mean decrease in accuracy mean decrease in MSE. Use names(datatrain) to see column numbers.
	imp <-as.data.frame(datatr.rf$importance[,1:2])
	wavenum <- rownames(imp) # get the waveband names
	wavenumd <- round(as.numeric(gsub("X","",wavenum)),0) # Get wavenumbers as numeric
	imp2 <- cbind(wavenumd,imp) # Add numeric wavenumbers
	
	
	## Plot Mean decrease accuracy %IncMSE against wavenumber. Deafult wavenmbers is 600,400.
	p1 <- ggplot(imp2,aes(imp2[,1],imp2[,2])) + geom_line() + scale_x_continuous(limits=c(500,4000)) + scale_x_reverse() + xlab(expression(paste("Wavenumber ","(",cm^-1,")"))) + ylab(expression(paste("%IncMSE"))) + 
	 opts(title = "", axis.text.x = theme_text(size=10),
	  axis.title.x = theme_text(hjust=0.55, vjust = 0, size=10),
	  axis.text.y = theme_text(size=10),
	  axis.title.y = theme_text(hjust=0.6, vjust=0.35,angle=90,size=10),
	  plot.title=theme_text(hjust=0.55,vjust=1,size=10,face="plain"),
	  panel.border=theme_rect(size=0.75),
	  legend.position="none"
	   )
	# Plot Mean decrease MSE IncNodePurity against wavenumber
	p2 <- ggplot(imp2,aes(imp2[,1],imp2[,3])) + geom_line() + scale_x_continuous(limits=c(500,4000)) + scale_x_reverse() + xlab(expression(paste("Wavenumber ","(",cm^-1,")"))) + ylab(expression(paste("IncNodePurity"))) + 
	 opts(title = "", axis.text.x = theme_text(size=10),
	  axis.title.x = theme_text(hjust=0.55, vjust = 0, size=10),
	  axis.text.y = theme_text(size=10),
	  axis.title.y = theme_text(hjust=0.6, vjust=0.35,angle=90,size=10),
	  plot.title=theme_text(hjust=0.55,vjust=1,size=10,face="plain"),
	  panel.border=theme_rect(size=0.75),
	  legend.position="none"
	   )
	# Print plots together
	grid.newpage()
	vport <- function(x, y) viewport(layout.pos.row=x,layout.pos.col=y);
	pushViewport(viewport(layout = grid.layout(nrow=2, ncol=1, widths=unit(6,"inches"), heights=unit(3,"inches"))))
	print(p1, vp=vport(1,1));
	print(p2, vp=vport(2,1));
	
	# Interactive plot of mean decrease MSE
	#spec.plot <- function(panel) 
	#{with(panel, {
	#with(imp2,{plot(imp2[,1],imp2[,2], type="l",xlim=rev(c(xmin,xmax)))}) 
	#minor.tick(nx=10)
	#})
		#panel
		#}
	#var.panel <- rp.control(xmin=600, xmax=4000,title="R panel")
	#rp.slider(var.panel,xmin,from=600,to=4000,showvalue=TRUE,action=spec.plot,title="x-axis Min")
	#rp.slider(var.panel,xmax,from=600,to=4000,showvalue=TRUE,action=spec.plot,title="x-axis Max")
	#############################################################################
	
	# Plot rel error vs number of trees
	plot(datatr.rf)
	
	# Feature selection using random cross-validation. Warning - takes some time!
	result <-rfcv(datatr[,-1], datatr[,1], cvfold=5)
	with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))
	
	# Frequency with which each variable is used. Use to find column no. of important variables
	varUsed(datatr.rf)
	   
	# Tune the mtry parameter
	tuneRF(x=datatr[,-1], y=datatr[,1], ntree=50)
	
	## Rank importance on %IncMSE (or node purity if choose column 3)
	imp3 <-imp2[order(imp2[,2], decreasing=TRUE),] # Sort importance variables in descending order 
	imp3 <-cbind(imp3,1:nrow(imp3)) # add row numbers as counter
	head(imp3,100)
	
	## Plot target variable value against absorbance for a given wavenmuber. Enter rank (nrank) of variable in imp3.
	nrank <- 1
	waveband <- paste(rownames(imp3[nrank,],)) # Get the name of the nth ranking waveband
	dataplot <- data.frame(datatr[yy],datatr[waveband]) # Subset y and x variable in new dataframe
	plottitle <- paste("Wavenumber", as.character(round(as.numeric(gsub("X","",rownames(imp3[nrank,]))),0)))
	p3 <- ggplot(dataplot, aes(dataplot[,2], dataplot[,1])) + geom_point() + xlab("Derivative absorbance") + ylab(yy) + 
	  opts(title = plottitle, axis.text.x = theme_text(size=11),
	  axis.title.x = theme_text(hjust=0.55, vjust = 0, size=12),
	  axis.text.y = theme_text(size=11),
	  axis.title.y = theme_text(hjust=0.6, vjust=0.35,angle=90,size=12),
	  plot.title=theme_text(hjust=0.55,vjust=1,size=13,face="bold"),
	  panel.border=theme_rect(size=0.75),
	  legend.position="none"
	   )
	 p3
	 
	 #### Plot top ranking wavenumbers as vertical lines on quartile mean raw spectrum #####################
	## Subset desired data and retain only complete cases
	#### Specify bands to highlight: 12 is main peak, 11 and 13 are margins
	b11 <- 0
	b12 <- 1687.4
	b13 <- 0
	b21 <- 0
	b22 <- 0
	b23 <- 0
	b31 <- 0
	b32 <- 0
	b33 <- 0
	## Prepare data
	data.raw <- data.frame(data.soil[,yy],data.spec) # subset the variables of interest
	data.raw2 <- data.raw[complete.cases(data.raw),] # retain only complete cases
	data.raw3 <- rename(data.raw2, c(data.soil...yy.="yyy")) # rename first col with dummy yyy
	qcut <- function(x, n) {
	  cut(x, quantile(x, seq(0, 1, length = n + 1)), labels = seq_len(n),
	    include.lowest = TRUE)
	}
	### Specify no of quartiles (nq). Ignore to four warning messages.
	nq<-3
	data.raw4 <- within(data.raw3, {quant <- qcut(data.raw3$yyy,eval(nq))})
	data.raw5 <- aggregate(data.raw4, by=list(data.raw4$quant), FUN=mean)
	data.sub <-  subset(data.raw5, select = c(X601.7:X4001.6)) 
	data.subt <- t(data.sub)
	data.subt1 <-as.data.frame(data.subt)
	data.subt1$wave <- rownames(data.subt1) # Get rownames as column
	data.subt1$wave <- round(as.numeric(gsub("X","",data.subt1$wave)),0) # Get wavenumbers as numeric; use raw waveband prefix
	colnames(data.subt1) <- paste("q",colnames(data.subt1),sep="") # Name the spectral columns
	data.subt2<- melt(data.subt1, id="qwave")
	
	## Plot data
	mspec <- ggplot(data.subt2, aes(x=qwave, y=value, colour=variable)) + geom_line() + xlab("wavenumber") + ylab("Absorbance") + scale_x_reverse(limits=c(4000,500)) +
	 geom_vline(xintercept = b11, colour="black", linetype = "dotted", size=0.35) + geom_vline(xintercept = b12, colour="black", linetype = "solid", size=0.35) + geom_vline(xintercept = b13, colour="black",  linetype = "dotted", size=0.35) +
	 geom_vline(xintercept = b21, colour="black",  linetype = "dotted", size=0.35) + geom_vline(xintercept = b22, colour="black",  linetype = "solid", size=0.35) + geom_vline(xintercept = b23, colour="black",  linetype = "dotted", size=0.35) +
	 geom_vline(xintercept = b31, colour="black",  linetype = "dotted", size=0.35) + geom_vline(xintercept = b32, colour="black",  linetype = "solid", size=0.35) + geom_vline(xintercept = b33,  colour="black",  linetype = "dotted", size=0.35) +
	  opts(title = "", axis.text.x = theme_text(size=11),
	  axis.title.x = theme_text(hjust=0.55, vjust = 0, size=12),
	  axis.text.y = theme_text(size=11),
	  axis.title.y = theme_text(hjust=0.6, vjust=0.35,angle=90,size=12),
	  plot.title=theme_text(hjust=0.55,vjust=1,size=13,face="bold"),
	  panel.border=theme_rect(size=0.75),
	  legend.position="right"
	   )
	mspec
	# Ignore one waveband beyond 4000 omitted
	
	 #### Plot top ranking wavenumbers as vertical lines on transformed spectrum #####################
	## Prepare data
	data.traw4 <- within(datatr, {quant <- qcut(datatr[,1],eval(nq))})
	data.traw5 <- aggregate(data.traw4, by=list(data.traw4$quant), FUN=mean)
	data.tsub <-  subset(data.traw5, select = c(X601.7:X4001.6)) 
	data.tsubt <- t(data.tsub)
	data.tsubt1 <-as.data.frame(data.tsubt)
	data.tsubt1$wave <- rownames(data.tsubt1) # Get rownames as column
	data.tsubt1$wave <- round(as.numeric(gsub("X","",data.tsubt1$wave)),0) # Get wavenumbers as numeric; use transformed waveband prefix
	colnames(data.tsubt1) <- paste("q",colnames(data.tsubt1),sep="") # Name the spectral columns
	data.tsubt2<- melt(data.tsubt1, id="qwave")
	#### Plot data
	mspect <- ggplot(data.tsubt2, aes(x=qwave, y=value, colour=variable)) + geom_line() + xlab("wavenumber") + ylab("Absorbance") + scale_x_reverse(limits=c(4000,500)) +
	 geom_vline(xintercept = b11, colour="black", linetype = "dotted", size=0.35) + geom_vline(xintercept = b12, colour="black", linetype = "solid", size=0.35) + geom_vline(xintercept = b13, colour="black",  linetype = "dotted", size=0.35) +
	 geom_vline(xintercept = b21, colour="black",  linetype = "dotted", size=0.35) + geom_vline(xintercept = b22, colour="black",  linetype = "solid", size=0.35) + geom_vline(xintercept = b23, colour="black",  linetype = "dotted", size=0.35) +
	 geom_vline(xintercept = b31, colour="black",  linetype = "dotted", size=0.35) + geom_vline(xintercept = b32, colour="black",  linetype = "solid", size=0.35) + geom_vline(xintercept = b33,  colour="black",  linetype = "dotted", size=0.35) +
	  opts(title = "", axis.text.x = theme_text(size=11),
	  axis.title.x = theme_text(hjust=0.55, vjust = 0, size=12),
	  axis.text.y = theme_text(size=11),
	  axis.title.y = theme_text(hjust=0.6, vjust=0.35,angle=90,size=12),
	  plot.title=theme_text(hjust=0.55,vjust=1,size=13,face="bold"),
	  panel.border=theme_rect(size=0.75),
	  legend.position="right"
	   )
	mspect 
	# Ignore one waveband beyond 4000 omitted
	
	#### Interactive plot of quartile spectra
	#spec.plot <- function(panel) {
	#with(panel, {
	#with(data.subt2,{plot(qwave,value, type="l",xlim=rev(c(xmin,xmax)))})
	#minor.tick(nx=5)
	#})
	#panel
	#}
	#var.panel <- rp.control(xmin=400, xmax=4000,title="R panel")
	#rp.slider(var.panel,xmin,from=400,to=4000,showvalue=TRUE,action=spec.plot,title="x-axis Min")
	#rp.slider(var.panel,xmax,from=400,to=4000,showvalue=TRUE,action=spec.plot,title="x-axis Max")
	#############################################################################
	
	# View a tree - not operational - only for classification?
	# treeset.randomForest(datatr.rf,5)
	
	########## Do not run beyond this point on first fit.
	
	## Rerun with important variables. Set number of trees. The re-run prediction plot again!
	nVar = 100 # Choose number of top ranking variables for model rerun
	newvar = rownames(imp3[1:nVar,]) # extract their waveband names
	datatr.rf <- randomForest(data=datatr, x=datatr[,newvar], y=datatr[,1], ntree=200, importance=TRUE, na.action=na.omit)
	datatr.rf
	
	setwd(mod.d)
	
	### Save the current model. Set the version no.
	version <- 1
	newmod = paste(yy,".rf_alpha",sep="")
	assign(newmod,datatr.rf)
	save(list=paste(newmod), file=paste(yy,".","v",version,".rfmodel",".Rdata",sep=""),compress="xz") # save with maxiumum compression
	
	## Print key results to files
	sink(paste(yy,".","v",version,".model",".txt",sep="")) # Start writing to file
	datatr.rf
	mp.train
	mp.oob
	imp3
	imp
	sink() # end writing to file
	
	## Write graph to png. (Don't adjust width & height as this changes font sizes).
	setwd(dia.d)
	png(file=paste(yy,".","v",version,".importance.png",sep=""),width=6,height=3, units = "in",res=300); 
	print(p1)
	dev.off()
	
	png(file=paste(yy,".","v",version,".response.png",sep=""),width=4,height=4, units = "in",res=300); 
	print(p3)
	dev.off()
	
	png(file=paste(yy,".","v",version,".bandplot.png",sep=""),width=5,height=4, units = "in",res=300); 
	print(mspec)
	dev.off()
	
	png(file=paste(yy,".","v",version,".bandplot_der.png",sep=""),width=5,height=4, units = "in",res=300); 
	print(mspect)
	dev.off()
	
	png(file=paste(yy,".","v",version,".predicted.png",sep=""),width=8,height=4, units = "in",res=300); 
	grid.newpage()
	vport <- function(x, y) viewport(layout.pos.row=x,layout.pos.col=y);
	pushViewport(viewport(layout = grid.layout(nrow=1, ncol=2, widths=unit(4,"inches"), heights=unit(4,"inches"))))
	print(mp1, vp=vport(1,1));
	print(mv1, vp=vport(1,2));
	dev.off()
}
######### Diagnostic test ######### DiagnosisMed library required 
rocdat <- data.frame(datatr[,yy], datatr.rf$predicted)
rocdat <- rename.vars(rocdat, c("datatr...yy.","datatr.rf.predicted"), c("gold","testval"))
hist(rocdat$gold) # Histogram y variable

#### Specify yclass cut-off; code case as 1 and non-case as 0
cutval <- 8

#### Use if values above cut off are case.
rocdat$goldbin <- ifelse(rocdat$gold > cutval, 1,0)
rocdat$testbin <- ifelse(rocdat$testval > cutval, 1,0)

#### Use if values below cut off are case
rocdat$goldbin <- ifelse(rocdat$gold < cutval, 1,0)
rocdat$testbin <- ifelse(rocdat$testval < cutval, 1,0)
rocdat$testval <- rocdat$testval*-1 # To correct the ROC direction

## Show diagnostic ROC curve and test results
rocplot <- with(rocdat,ROC(goldbin,testval,Plot.point="Error.rate" )) # Standard plot
diagnos <- with(rocdat,diagnosis(goldbin,testbin)) # Test results
# Interactive plot
#with(rocdat,interact.ROC(goldbin,testval)) # Interactive plot

rocres <- rocplot$test.diag.table
rocres.min <-  rocres[order(rocres$Error.rate),]
min.error.stats <- rocres.min[1,]
min.error.stats
rocres.max <- rocres[order(rocres$Efficiency, decreasing=TRUE ),]
max.eff.stats <- rocres.max[1,]
max.eff.stats

# Write test results to file
setwd(mod.d)
sink(paste(yy,".diagnostic.",cutval,".txt",sep="")) 
diagnos
min.error.stats
max.eff.stats
sink()

#############################
