require(caret)
require(soil.spec)
require(readr)

	#STEP1.Read new spectra with unknown soil property
	new_spectra<-read_csv("")#This is a fast way of reading large spectral files
		#Note1: Ensure the new table contains the same wavebands as the spectra table used for developing calibration models
		#Note2: Preprocess the new spectra using sgolay first derivative method if not already not done so.

	#STEP2. Read the model to be used e.g for clay
	mir.soil.property1<-readRDS("")#add full path to your model e.g for  my clay model in d:/ Models/Clay.RDS the full path is  "d:/ Models/Clay.RDS
		#Note 3: If you have n models read each at a time and name them appropriately. e.g if I have clay and carbon models I will name them as:
		
	#model1: clay<-readRDS(" d:/ Models/Clay.RDS")
	#model2: carbon<-readRDS(" d:/ Models/Clay.RDS")
	
	#STEP3: using the new_spectra after processing,predict each soil property from the same.
	predicted_clay<-predict(clay,new_spectra)
	predicted_carbon<-predict(carbon, new_spectra)
	
	#Combine the predicted soil properties to get new_spectra_id, predicted_clay & predicted_carbon as:
	new_spectra_id<-new_spectra[,k]#k is the column with spectra id normaly column 1.
	predicted<-cbind(new_spectra_id, predicted_clay, predicted_carbon)
	
	#Write predicted into a csv/txt file
	write.table(predicted,"MIR predicted soil properties.csv", sep=",",row.names=TRUE)
