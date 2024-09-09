######################################################################
# carv.R: R functions for generating Circadian Activity Rhythm Vars
# (c) Rosemary Braun <rbraun@northwestern.edu> 2017
######################################################################
#
# This script is a conversion of the SAS code to generate 
# Circadian Activity Rhythm Variables.
#
# For background, see files:
# - "Instructions on how to generate Circadian Activity Rhythm Variables (1).pdf"
# included in mail from Kathy Reid on Wed 28 Jun 2017 (Subject: RE: Timestamp help).
#
# Usage:
#
# - Place all the actiware .csv files you wish to analyze in a single
#   folder.  You need not strip the headers or manipulate them in any
#   way.  Be sure nothing else is in that folder, and make note of 
#   the folder location; this will be your PATH. (quote it)
#
# - Source this R script:
#     source("carv.R")
#
# - Use the wrapper functions to read CSVs & run CARV analysis:
#     my.CSV.files <- dir(PATH, full.names=T)
#     actiware.data <- readActiwareCSVs(my.CSV.files)
#     my.linearFit <- lActReg(actiware.data)
#     my.nonlinearFit <- fitNLS(actiware.data)


#     NLSredux(actiware.data)
#     # plot with
#     plotfit2(ID,actiware.data,my.nonlinearFit)
#
# - In the above,
#     path = character string giving the path to the actiware files
#     ID   = character string for a particular subject/file for whom
#            you wish to plot the nonlinear fit  

#=====================================================================
# Ancillary functions
#=====================================================================

.timeString24 <- function(timeStrings){
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


#=====================================================================
# Reading in the data
#=====================================================================

readActiwareCSV <- function(filename,idcode=filename){
  tmp <- readLines(filename)
  
  #**Find the subject**************************************************
  idRow <- grep("Identity:",tmp)
  
  Id_data1<- read.csv(filename, skip=  idRow-1, nrows = 1, header = F)
  
  Id_data<-Id_data1[2]
  
  ##******************************************************************
  
  
  headRow <- grep("Line.*Activity",tmp)
  if(length(headRow)>1){
    stop("Input file misformatted; >1 section with Activity data")
  }
  cat("- Reading", filename,"... \n")
  if(length(headRow)==0){
    # We have read in a file with the stripped header;
    # this code follows that of the SAS program
    out <- read.csv(filename,header=FALSE,stringsAsFactors=F)
    colnames(out) <- c("line","dateString","timeString","Activity","light","slpwk","interval")
    out$time <- .timeString24(out$timeString)
    out$lumeday <- as.Date(out$dateString,"%m/%d/%Y")
    # computed in the SAS code but never used:
    #   out$outbed <- 0
    #   out$outbed[toupper(out$interval)=="ACTIVE"] <- 1
    #   out$inbed <- 0
    #   out$inbed[toupper(out$interval)=="REST"] <- 1
    #   out$inbed[toupper(out$interval)=="REST-S"] <- 1
    #   out$autoact <- 1-out$slpwk
    # set up a clock time variable that counts starting on the first lumeday
  } else {
    # We have a raw Actiware .csv file;
    # read it in
    out <- read.csv(filename,skip=(headRow-1),quote="\"",stringsAsFactors=F) 
    
    if(out$Time[1] ==""){out<-out[-1,]} # to avoide first row for some formated csv file
    
    out$dateString <- out$Date
    out$timeString <- out$Time
  }
  
  ui_ind<-ind_st.time
  idRow2 <- grep(ui_ind, out$timeString)
  out<-out[c(idRow2[1]:nrow(out)), ]
  
  
  out$time <- .timeString24(out$timeString)
  out$lumeday <- as.Date(out$dateString,"%m/%d/%Y")
  first.lumeday <- min(as.numeric(out$lumeday))
  out$cloktime <- 24*(as.numeric(out$lumeday)-first.lumeday)+out$time
  # Add idcode and phs...
  # "phs" was in the SAS code as an indicator of study phase, 
  # but always kept at 0; should probably be removed in future versions
  
  out$phs <- 0
  out$idcode <- as.character(idcode)
  out$sub.ID <- as.character(Id_data$V2)# Subject ID
  
  
  
  
  # chrono sort the output
  
  out <- out[order(out$cloktime),]
  rownames(out) <- NULL
  out$Activity[out$Activity<0] <- NA
  
  

  
  
  out$lActivity <- log10(out$Activity+1)
  
  return(out[,c("idcode","sub.ID", "phs","time","Activity","cloktime","lActivity")])
}






readActiwareCSVs <- function(filenames,idcodes=basename(filenames)){
  if(length(filenames)!=length(idcodes)){
    stop("Need an idcode for each file; filenames and idcodes should have same length")
  }
  cat("Reading",length(filenames),"files...\n")
  out <- lapply(1:length(filenames),function(i){readActiwareCSV(filenames[i],idcodes[i])})
  out <- do.call(rbind,out)
  return(out)
}

# example:
# tmp <- readActiwareCSVs(dir("testData",pattern="00",full.names=T),c("0008","1008","2008"))


#=====================================================================
# Linear model
#=====================================================================


lActReg <- function(actidat){
  subjcode <- paste(actidat$idcode,actidat$phs,sep=".")
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
    try(nls(
      lActivity ~ curveFn(cloktime,actamp,actbeta,actphi,actmin,actalph),
      data = actidat[actidat$idcode==ID,],
      start = starts, lower = lowerB, upper = upperB, algorithm="port"
    ))
  })
  names(out) <- as.character(linfit$idcode)
  return(out)
}







NLSredux <- function(actidat){
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
  
  
  out$missing.data <- sapply(split(actidat$lActivity,actidat$idcode),function(z){
    a_miss<-aggr(z)
    round(a_miss$percent[2], digits = 2)
  }) 
  
  
  

  
  
  out <- out[,c("idcode","missing.data","actamp","actbeta","actphi","actmin","actmesor","actupmesor","actdownmesor","actalph","actwidthratio","RsqNLS","Fact","Fnlrgact")]
  return(out)
}






#=====================================================================
# Smoothing; cf DARPA\ Timestamp\ Actigraphy\ Scoring\ 4.10.17.docx
#=====================================================================
smoothAct30 <- function(actidat){
  # assumes 30-sec epochs
  smoothAct30sub <- function(actisub,sepoch=30){
    actisub <- actisub[order(actisub$cloktime),]
    epoch.no <- actisub$cloktime*60*60/sepoch
    epoch.no <- epoch.no-min(epoch.no)+1+4 # buffer of 4 behind for 30s epochs 
    epoch.last <- max(epoch.no)+4 # buffer of 4 ahead for 30s epochs
    epoch.id <- paste("epoch", epoch.no,sep=".")
    epoch.act <- actisub$Activity
    names(epoch.act) <- epoch.id
    smooth.act <- rep(NA,epoch.last) # buffer of 4 ahead & 4 behind for 30s epochs
    names(smooth.act) <- paste("epoch", 1:epoch.last,sep=".")
    smooth.act[names(epoch.act)] <- epoch.act
    shift.act <- rbind(
      smooth.act[(1:(length(smooth.act)-8)+0)],
      smooth.act[(1:(length(smooth.act)-8)+1)],
      smooth.act[(1:(length(smooth.act)-8)+2)],
      smooth.act[(1:(length(smooth.act)-8)+3)],
      smooth.act[(1:(length(smooth.act)-8)+4)],
      smooth.act[(1:(length(smooth.act)-8)+5)],
      smooth.act[(1:(length(smooth.act)-8)+6)],
      smooth.act[(1:(length(smooth.act)-8)+7)],
      smooth.act[(1:(length(smooth.act)-8)+8)]
    )
    colnames(shift.act) <- names(smooth.act[(1:(length(smooth.act)-8)+4)])
    smooth.coef <- c(1/25,1/25,1/5,1/5,2,1/5,1/5,1/25,1/25)
    smooth.num <- colSums(smooth.coef*shift.act,na.rm=T)
    # normalize? apparently not per manual
    smooth.den <- colSums(smooth.coef*!is.na(shift.act),na.rm=T)
    actisub$smoothact <- smooth.num[names(epoch.act)]
    actisub$lsmoothact <- log10(actisub$smoothact+1)
    return(actisub)
  }
  out <- do.call(rbind(lapply(split(actidat,actidat$idcode), smoothAct30sub)))
  return(out)
}
  
  
#=====================================================================
# plotting routines... needs cleanup!
#=====================================================================

plotfit <- function(id,actidat,nlsfit){
  #args <- as.list(summary(nlsfit[[id]])$coef[,1])
  #args$cloktime <- actidat[actidat$idcode==id]$cloktime
  plotDF <- actidat[actidat$idcode==id,c("cloktime","lActivity")]
  plotDF <- na.omit(plotDF)
  plotDF$predlact <- predict(nlsfit[[id]],plotDF$cloktime)
  ptcol <- c("black","grey")[(plotDF$cloktime%/%24)%%2+1]
  plot(plotDF$cloktime, plotDF$lActivity,main=id,cex=0.7,col=ptcol,xlab="cloktime",ylab="log10(Activity+1)")
  lines(plotDF$cloktime, plotDF$predlact,col=2,lwd=2)
}


plotfit2 <- function(id,actidat,nlsfit,type="p",daysPerPage=4){
  #args <- as.list(summary(nlsfit[[id]])$coef[,1])
  #args$cloktime <- actidat[actidat$idcode==id]$cloktime
  plotDF <- actidat[actidat$idcode==id,c("time","cloktime","lActivity")]
  plotDF <- na.omit(plotDF)
  plotDF$day <- (plotDF$cloktime%/%24)+1
  plotDF$predlact <- predict(nlsfit[[id]],plotDF$cloktime)
  plotDF <- split(plotDF,plotDF$day)
  old.par <- par(no.readonly=TRUE)
  on.exit({par(old.par);par(mfrow=c(1,1))})
  par(mfrow=c(daysPerPage,1),mar=c(3,5,1,1))
  if(interactive()&length(plotDF)>daysPerPage){
    par(ask=T)
  }
  for (day in 1:length(plotDF)){
    plot(plotDF[[day]]$time, plotDF[[day]]$lActivity, xlim=c(0,24),
      main=paste(id,"- day",day), xlab="cloktime",ylab="log10(Activity+1)",type="n")
    polygon(c(0,6,6,0),c(0,0,10000,10000),col="grey",border=0)
    polygon(c(20,24,24,20),c(0,0,10000,10000),col="grey",border=0)
    lines(plotDF[[day]]$time, plotDF[[day]]$lActivity,type=type,cex=0.5)
    lines(plotDF[[day]]$time, plotDF[[day]]$predlact,col=2,lwd=2)
  }
}


#=====================================================================
# Example usage -- needs documentation -- also see "Usage" at top
#=====================================================================

# if(FALSE){ # example:
#   rm(list=ls())
#   source("carv.R")
#   tmp <- readActiwareCSVs(dir("testData",pattern="00",full.names=T),c("0008","1008","2008"))
#   #tmp2 <- readActiwareCSVs(dir("/Volumes/fsmresfiles/Neurology/CRSR/Projects/Timestamp - DARPA/Actigraphy/QCed CSV FILES",pattern="1741",full.names=T))
#   foo <- lActReg(tmp)
#   foo2 <- fitNLS(tmp)
#   NLSredux(tmp)
#   plotfit2("2008",tmp,foo2)
# }
