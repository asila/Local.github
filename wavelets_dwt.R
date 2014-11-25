require(wavelets)
##Read file with raw spectra for wavelet transformation
raws<-read.csv(file.choose())#The table read should contain spectra in all the columns except in column
									#column 1 contains spectra ids 
		raw<-raws[,-1]				#remove first column
		waveb <- as.numeric(substr(colnames(raw),2,19))
        waveb.1024.up <- round(max(waveb))
        waveb.1024.down <- round(min(waveb))
        n<-as.character(c(31,127,255,511,1023))
        waven<-menu(n,graphics=TRUE,title="Select columns to make spectra to power 2")
        waveb.1024.n <-as.numeric(n[waven])
        waveb.1024.step <- (waveb.1024.up - waveb.1024.down)/waveb.1024.n
        waveb.1024 <- c()
        waveb.1024[1] <- waveb.1024.down
        for (i in 2: (waveb.1024.n+1)) {
            waveb.1024[i] <- round(waveb.1024.down + ((i - 1) * waveb.1024.step), 5)
        }
                rownames(raw)<-raws[,1]
        raw.comp <- matrix(nrow = nrow(raw), ncol = length(waveb.1024), 
            dimnames = list(rownames(raw), waveb.1024))
        for (i in 1:nrow(raw)) {
            raw.comp[i, ] <- round(spline( waveb, raw[i, ], method = "natural", 
                xout = waveb.1024)[[2]], 6)
        }
        lev <- (waveb.1024.n+1)/4
        slo <- 3
        filte = "haar"
        transdf <- matrix(nrow = nrow(raw.comp), ncol = lev, 
            dimnames = list(rownames(raw.comp), paste("WC_", 
                c(1:lev), sep = "")))
        for (i in 1:nrow(transdf)) {
            blub <- dwt(raw.comp[i, ], filter = filte)
            transdf[i, ] <- slot(blub, "W")[[slo]]
        }
transdf<-cbind(as.vector(row.names(transdf)),transdf)
colnames(transdf)<-c("SSN",colnames(transdf[,-1]))
write.table(transdf,file="Wavelets5.csv",sep=",",row.names=FALSE)