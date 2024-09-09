



stopifnot(
  require(nparACT)
)

## This is is the program for running carv and nparact by merging different files


.timeString24 <- function(timeStrings){
  processSplit <- function(timeSplit){
    hr <- as.numeric(timeSplit[1])
    min <- as.numeric(timeSplit[2])
    sec <- as.numeric(timeSplit[3])
    isPM <- toupper(substr(timeSplit[4],1,1))=="P" # check whether is it PM ar Am
    if(hr==12){isPM <- !isPM}
    hr <- (hr+12*isPM)%%24
    return(hr+min/60+sec/3600) # Convert into 24hr format ....... Convert 10:30:30 ------> 10.5083
  }
  
  
  tryCatch(timeSplit <- strsplit(timeStrings,"[: ]"), # split "12:00:00 PM" -------> "12" "00" "00" "PM",
           
           error = function(e) {
             stop('Subject not found')})
  
  
  
  
  sapply(timeSplit, processSplit)
}

## *******************************************************************************************



.nparACT.timeString24 <- function(timeStrings){
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




#************************************************************************************************


########### NEW FUNCTION ##############


#*********************************************************************************************

#readActiwareCSV ()
# Remove all unwanted item from csv file, provide header and data only
# This function calls .timeString24


#*********************************************************************************************


readActiwareCSV <- function(filename,idcode=filename){
  
  out <- filename
  out$time <- .timeString24(out$timeString) # convert time into 24 hr format 10:30:30 ------> 10.5083
  out$lumeday <- as.Date(out$dateString,"%m/%d/%Y")
  first.lumeday <- min(as.numeric(out$lumeday)) # get the first day
  out$cloktime <- 24*(as.numeric(out$lumeday)-first.lumeday)+out$time # continues timing 24, 24.1,....48...72.....
  
  
  #Day # time     clock time
  # 1    12.00   00.00
  #1    23.9    23.9
  #2    0       24
  #2    1       25
  
  
  out$phs <- 0
  
  out$idcode <- as.character("AB")
  
  

  
  # ***************************************************
  out$npartime <- .nparACT.timeString24(out$timeString)
  out$DayTime <- paste(out$lumeday, out$npartime)
  # ***************************************************
  
  
  # chrono sort the output
  
  out <- out[order(out$cloktime),] # sort increasing order of $clocktime
  rownames(out) <- NULL
  out$Activity[out$Activity<0] <- NA
 
  
  ##************************************************
  start.date1=as.Date(st.day,"%m/%d/%Y")
  end.date1=as.Date( ed.day,"%m/%d/%Y")
  D1=as.numeric(start.date1)
  D2=as.numeric(end.date1)




  T1=.timeString24(st.time)
  T2=.timeString24(ed.time)
  
  out$cloktime <- 24*(as.numeric(out$lumeday)-first.lumeday)+out$time 
  TC1 <- 24*(D1-first.lumeday)+T1
  TC2 <- 24*(D2-first.lumeday)+T2
  ##************************************************  
  
  #out<-subset(out,  as.numeric(out$lumeday)>=D1 & as.numeric(out$lumeday)<=D2)
  
  out<-subset(out,    out$cloktime >=TC1 & out$cloktime <= TC2)
  
  
  ##*********************************************************************************************
 # Down sampling 
  
  T_samp=as.numeric(Ts)/as.numeric(default.sampling.time1) # Find how many samples has to skip
  
  x1=out$Activity
  x2=out$White.Light
  
  
  x1m=matrix(x1, nrow = T_samp) # split the data
  x2m=matrix(x2, nrow = T_samp)
  
  x1na=apply(x1m, 2, pct_na) # Find fraction of NaN in each splited group
  x2na=apply(x2m, 2, pct_na)
  
  Ix1=which(x1na>=0.5) # find the index 
  Ix2=which(x2na>=0.5)
  
  DS_data1<- colMeans(matrix(x1, nrow = T_samp), na.rm=TRUE) # Mean of each splited group
  DS_data2<- colMeans(matrix(x2, nrow = T_samp), na.rm=TRUE)
  
  DS_data1[Ix1]<-NA  # replace mean with NaN if toatal NaN > 50%
  DS_data2[Ix2]<-NA

  
  out<-out[seq(1,nrow(out),T_samp), ] # Down sampling (Skipping T_sample)
  
  out$Activity<-DS_data1[1:nrow(out)]  # replace activity with average
  out$White.Light<-DS_data2[1:nrow(out)]
  
  out$lActivity <- log10(out$Activity+1) # find the log of (Activit+1)  ------>  Activit=12, lActivity=log(13)
  
  ##*************************************************************************************************
  
  
  
  
  return(out[,c("idcode", "phs","time","DayTime", "Activity","cloktime","lActivity", "White.Light", "Time", "dateString")]) # output of readActiwareCSV function
  
  
}


#Eg-:readActiwareCSV()



Merge.ActiwareCSV <- function(filename){
  
 
  #modified on 5/10/2020
  
  
  tmp <- readLines(filename)
  
  #**Find the subject**************************************************
  idRow <- grep("Identity:",tmp)
  
  Id_data1<- read.csv(filename, skip=  idRow-1, nrows = 1, header = F)
  
  Id_data<-Id_data1[2]
  
  ##******************************************************************

  
  headRow <- grep("Line.*Activity",tmp) # Give the row nmber where "Activity" header is present
  
  if(length(headRow)>1){
    stop("Input file misformatted; >1 section with Activity data")
  }
  
  
  cat("- Reading", filename,"... \n")
  
  
  if(length(headRow)==0) 
  {
    # We have read in a file with the stripped header;
    # this code follows that of the SAS program
    
    out <- read.csv(filename,header=FALSE,stringsAsFactors=F)
    colnames(out) <- c("line","dateString","timeString","Activity","light","slpwk","interval") # set the column name as.....
    out$time <- .timeString24(out$timeString)
    out$lumeday <- as.Date(out$dateString,"%m/%d/%Y")
    
  } else {
    
    
    # We have a raw Actiware .csv file;
    # read it in
    
    
    out <- read.csv(filename,skip=(headRow-1),quote="\"",stringsAsFactors=F) # Remove the line above the true headings (Now only data)
   
    
    if(out$Time[1] ==""){out<-out[-1,]} # to avoide first row for some formated csv file
    
     out$dateString <- out$Date
    out$timeString <- out$Time
    out$sub.ID <- as.character(Id_data$V2)# Subject ID
    
    return(out[,c("sub.ID", "Date", "Time", "Activity", "White.Light", "dateString", "timeString")]) # output of readActiwareCSV function
    
  
  }
  

}


#Eg:-Merge.ActiwareCSV(dir(PATH, full.names=T))




           ##################################################################################################

                                           #lActReg ()

          #*******************************************************************************************


lActReg <- function(actidat){
  subjcode <- paste(actidat$idcode,actidat$phs,sep=".") # Concatenate vectors after converting to character
  # subjectcode= chr[1:length(actidata)], each element is "filename.csv.0" 
  
  actilist <- split(actidat,subjcode)
  
  fit.linact <- function(actisubj){
    actisubj$xcos <- cos((2*pi/24)*actisubj$cloktime)
    actisubj$xsin <- sin((2*pi/24)*actisubj$cloktime)
    regfit <- (lm(lActivity~xcos+xsin,data=actisubj))
    regout <- summary(regfit)
    out <- data.frame(
      idcode = actisubj$idcode[1],
      phs = actisubj$phs[1],
      linact.adjR2 = regout$adj.r.squared,
      linactxb = regout$coef[1,1],
      linactcos = regout$coef[2,1],
      linactsin = regout$coef[3,1],
      linactRSS = sum(regfit$residuals^2,na.rm=T),
      stringsAsFactors = FALSE
    )
    out$linactamp <- sqrt(out$linactcos^2+out$linactsin^2)
    out$linactmin <- out$linactxb-out$linactamp
    phase <- (atan2(out$linactsin,out$linactcos)+2*pi)%%(2*pi)
    out$linactacro <- phase*24/(2*pi)
    
    
    return(out)
  }
  actres1 <- do.call(rbind, lapply(actilist, fit.linact))
  return(actres1)
}




#=====================================================================
# The curve...
#=====================================================================

curveFn <- function(
  cloktime=seq(0,48,len=48*60+1), 
  actamp=0.001, # 0 < actamp < 5
  actbeta=2.00, # 0 < actbeta
  actphi=12,    # -3 < actphi < 27
  actmin=0,     # 0 <= actmin
  actalph=0.0   # -1 <= actalph <= 1 
){
  pi12 <- 2*pi/24
  # following the SAS code conventions...
  rhythm <- cos(pi12*(cloktime-actphi))
  expt <- exp( actbeta*(rhythm-actalph) )
  er <- expt/(1+expt)
  xmodel <- actmin+actamp*er
  # zmodel <- actmin+actamp*(exp(actbeta*(cos(pi12*(cloktime-actphi))-actalph))/(1+exp(actbeta*(cos(pi12*(cloktime-actphi))-actalph))))
  return(xmodel)
}


#=====================================================================
# The NLS fit...
#=====================================================================


fitNLS <- function(actidat){
  linfit <- lActReg(actidat)
  #actidat <- split(actidat,idcode)
  #linfit <- split(linfit,idcode)
  out <- lapply(linfit$idcode, function(ID){
    starts <- c(
      # actamp=0.001, # 0 < actamp < 5
      actamp=5*linfit$linactamp[linfit$idcode==ID],
      actbeta=2.00, # 0 < actbeta
      # actphi=12,    # -3 < actphi < 27
      actphi=linfit$linactacro[linfit$idcode==ID],
      # actmin=0, # 0 <= actmin
      actmin=max(linfit$actmin[linfit$idcode==ID],0),
      actalph=0.0   # -1 <= actalph <= 1 
    )
    lowerB <- c(actamp=0,actbeta=0,actphi=-3,actmin=0,actalph=-1)
    upperB <- c(actamp=5,actbeta=Inf,actphi=27,actmin=Inf,actalph=1)
    starts <- pmax(lowerB,starts)
    starts <- pmin(upperB,starts)
    starts <- as.list(starts)
    
    try(curve.model <-nls(
      lActivity ~ curveFn(cloktime,actamp,actbeta,actphi,actmin,actalph),
      data = actidat[actidat$idcode==ID,],
      start = starts, lower = lowerB, upper = upperB, algorithm="port"
    ))
    
  })
  
  names(out) <- as.character(linfit$idcode)
  return(out)
}




#********************************************************************************************

## NLSredux() 

# ********************************************************************************************



NLSredux <- function(actidat){
  
  # Finding the missing data
  

  linfit <- lActReg(actidat)  # this gets executed 2x -- very inefficient!
  nlsfit <- fitNLS(actidat)
  out <- as.data.frame(t(sapply(nlsfit,coef)))
  out$idcode <- rownames(out)
  actacos <- acos(out$actalph)/(pi/12)
  actacos[abs(out$actalph)>0.99] <- 6
  out$actupmesor <- -actacos+out$actphi
  out$actdownmesor <- actacos+out$actphi
  out$actwidthratio <- actacos/12
  out$n <- sapply(nlsfit,function(z){sum(summary(z)$df,na.rm=T)})
  out$RSS <- sapply(nlsfit,function(z){sum(residuals(z)^2,na.rm=T)})
  out$totSS <- sapply(split(actidat$lActivity,actidat$idcode),function(z){
    sum((z-mean(z,na.rm=T))^2,na.rm=T)
  }) 
  out$RsqNLS <- (out$totSS-out$RSS)/out$totSS
  out$Fact <- ((out$totSS-out$RSS)/4) / (out$RSS /(out$n-5))
  out$Fnlrgact <- ((linfit$linactRSS-out$RSS)/2) / (out$RSS /(out$n-5))
  out$actmesor <- out$actmin+out$actamp/2
  #out <- out[,c("idcode","actamp","actbeta","actphi","actmin","actmesor","actupmesor","actdownmesor","actalph","actwidthratio","RsqNLS","Fact","Fnlrgact")]
  
  a_miss <-  aggr(actidat$Activity)
  out$missing.data<-round(a_miss$percent[2], digits = 2)
  
  out <- out[,c("missing.data", "actamp","actbeta","actphi","actmin","actalph", "actmesor","actupmesor","actdownmesor","actwidthratio","RsqNLS","Fact","Fnlrgact")]
  #out <- out[,c("actamp","actbeta","actphi","actmin","actalph", "actmesor","actupmesor","actdownmesor","actwidthratio","RsqNLS","Fact","Fnlrgact")]
  return(out)
}






































