######################################################################
# nparactiware.R: 
#   wrapper to apply the nparact to actiware output files
# (c) Rosemary Braun <rbraun@northwestern.edu> 2018
######################################################################
#
# This file defines two wrapper functions for the nparACT package:
#
# * mod_nparACT_base_loop reads in actiware .csv files and runs
#   the nparACT_base_loop function from nparACT;
#
# * mod_nparACT_flex_loop reads in actiware .csv files and runs
#   the nparACT_flex_loop function from nparACT.
#
# Additionally, the "rewriteActiwareCSVs" takes in actiware .csv
# files and writes them out in the format nparACT expects for 
# its input.  YOU NEED NOT RUN THIS BY HAND (it is automatically
# called from the mod_nparACT_*_loop functions), but it is available
# should you wish to do so.

# Usage:
#
# - Place all the actiware .csv files you wish to analyze in a single
#   folder.  You need not strip the headers or manipulate them in any
#   way.  Be sure nothing else is in that folder, and make note of 
#   the folder location; this will be your path.
#
# - Source this R script:
#     source("nparactiware.R")
#
# - Use the wrapper functions to read CSVs & run nparACT analysis:
#     mod_nparACT_base_loop(path, outDir, SR, cutoff, plot, fulldays)
#     mod_nparACT_flex_loop(path, outDir, SR, cutoff, minutes, plot, fulldays)
#
# - In the above,
#     path   = character string giving the path to the actiware files
#     outDir = folder where you want the "rewritten for nparACT" 
#              files to go; by default, a temp directory. 
#     SR, cutoff, minutes, plot, fulldays, ... 
#            = these are all arguments to the nparACT functions; see
#              help pages for nparACT_base_loop & nparACT_flex_loop


#=====================================================================
# Package dependencies
#=====================================================================

stopifnot(
	require(nparACT)
)

#=====================================================================
# Ancillary functions
#=====================================================================

.timeDec24 <- function(timeStrings){
	processSplit <- function(timeSplit){
		hr <- as.numeric(timeSplit[1])
		min <- as.numeric(timeSplit[2])
		sec <- as.numeric(timeSplit[3])
		isPM <- toupper(substr(timeSplit[4],1,1))=="P"
		if(hr==12){isPM <- !isPM}
		hr <- (hr+12*isPM)%%24
		return(hr+min/60+sec/3600)
	}
	timeSplit <- strsplit(timeStrings,"[: ]")
	sapply(timeSplit, processSplit)
}


.timeString24 <- function(timeStrings){
	processSplit <- function(timeSplit){
		hr <- as.numeric(timeSplit[1])
		min <- as.numeric(timeSplit[2])
		sec <- as.numeric(timeSplit[3])
		isPM <- toupper(substr(timeSplit[4],1,1))=="P"
		if(hr==12){isPM <- !isPM}
		hr <- (hr+12*isPM)%%24
		timestr <- sprintf("%02.0f:%02.0f:%02.0f",hr,min,sec)
		return(timestr)
	}
	timeSplit <- strsplit(timeStrings,"[: ]")
	sapply(timeSplit, processSplit)
}

#=====================================================================
# Reading in the data
#=====================================================================

.ingestActiwareCSV <- function(filename,idcode=filename){
	tmp <- readLines(filename)
	headRow <- grep("Line.*Activity",tmp)
	if(length(headRow)>1){
		stop("Input file misformatted; >1 section with Activity data")
	}
	cat("- Reading", filename,"... \n")
	if(length(headRow)==0){ # stripped header
		out <- read.csv(filename,header=FALSE,stringsAsFactors=F)
		# NOTE HARD-CODED ASSUMPTION REGARDING COLUMNS for stripped headers
		colnames(out) <- c("line","dateString","timeString","Activity","light","slpwk","interval")
	} else { # raw Actiware .csv file with headers
		out <- read.csv(filename,skip=(headRow-1),quote="\"",stringsAsFactors=F) 
		
		if(out$Time[1] ==""){out<-out[-1,]} # to avoide first row for some formated csv file
		
		
		out$dateString <- out$Date  # rename col
		out$timeString <- out$Time  # rename col
	}
	

	
	ui_ind<-ind_st.time
	idRow2 <- grep(ui_ind, out$timeString)
	out<-out[c(idRow2[1]:nrow(out)), ]
	

	
	out$time <- .timeString24(out$timeString)
	out$dectime <- .timeDec24(out$timeString)
	out$day <- as.Date(out$dateString,"%m/%d/%Y")
  out$DayTime <- paste(out$day, out$time)
	first.day <- min(as.numeric(out$day))
	out$abstime <- 24*(as.numeric(out$day)-first.day)+out$dectime
	# Add idcode and phs...
	out$idcode <- as.character(idcode)
	# chrono sort the output
	out <- out[order(out$abstime),]
	rownames(out) <- NULL
	out$Activity[out$Activity<0] <- NA
	
	
	


	out$lActivity <- log10(out$Activity+1)
	return(out[,c("idcode","DayTime","time","Activity","abstime","lActivity")])
}






rewriteActiwareCSVs <- function(path,idcodes=path,outDir=file.path(tempdir(),"nparactiware")){
	#filenames <- dir(path,full.names=T)
  
  
  
	filenames <- path
	
	idcodes=seq(1:length(filenames))
	dir.create(outDir)
	# if(length(filenames)!=length(idcodes)){
	# 	stop("Need an idcode for each file; filenames and idcodes should have same length")
	# }
	cat("Reading",length(filenames),"files...\n")
	
	
	
	for(i in 1:length(filenames)){
	  
		tmp <- .ingestActiwareCSV(filenames[i],idcodes[i])[,c("DayTime","Activity")]
		tmp <- na.omit(tmp)
	#	write.table(tmp,file=file.path(outDir,paste("PROCESSED-",idcodes[i],sep="")),
		            write.table(tmp,file=file.path(outDir,paste("PROCESSED-",idcodes[i],sep="")),
			sep=",",row.names=F,col.names=F,quote=F)
		
	}
	
	
  cat("Files left in",outDir,"\n")
	return(outDir)
}



# example:
#  
# tmpout <- rewriteActiwareCSVs("testData2")

#=====================================================================
# Wrappers for nparACT loop functions
#=====================================================================

mod_nparACT_base_loop <- function(path, outDir=file.path(tempdir(),"nparactiware"), SR, cutoff = 1, plot = T, fulldays = T){
	processedFilePath <- rewriteActiwareCSVs(path, outDir=outDir)
	out <- nparACT_base_loop(processedFilePath, SR=SR, cutoff=cutoff, plot=plot, fulldays=fulldays)
	return(out)
}

mod_nparACT_flex_loop <- function(path, outDir=file.path(tempdir(),"nparactiware"), SR, cutoff = 1, minutes, plot = T, fulldays = T){
	processedFilePath <- rewriteActiwareCSVs(path, outDir=outDir)
	out <- nparACT_flex_loop(processedFilePath, SR=SR, cutoff=cutoff, minutes=minutes, plot=plot, fulldays=fulldays)
	return(out)
}

# example:
#
# mod_nparACT_base_loop("testData2", SR=2/60)
