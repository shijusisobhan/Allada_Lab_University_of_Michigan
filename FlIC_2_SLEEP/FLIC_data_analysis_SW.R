
# This is the all user function for FLIC data analysis


packages = c("ggplot2", "stats", "gridExtra", "reshape2","gtools")

## Load or install the required packges
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)




# require(ggplot2)
# require(stats)
# require(gridExtra)
# require(reshape2)
# require(gtools)


# ***************PriviteFunction*************************************************************************************

##### Basic calculation functions #####
## This function takes a vector of dates (as strings), a vector
## of times (24 hour as string) and a parameters object.
## it returns the elapsed seconds.
GetElapsedSeconds<-function(dfm){
  dates<-dfm$Date
  times<-dfm$Time
  ms <-dfm$MSec
  ## For US culture.
  tmp<-as.character(dfm$Time[1])
  tmp2<-regexpr('.M',tmp)
  if(tmp2[1]>-1) {
    ## There seems to be some confusion about the nature of the time stamp
    ## from different MCU. 
    ## Use this one if time in AM/PM
    fulltimes<-as.POSIXct(paste(dates,times),format="%m/%d/%Y %I:%M:%S %p",tz='UTC')
  }
  else {
    ## Use this one if time is military time.
    fulltimes<-as.POSIXct(paste(dates,times),format="%m/%d/%Y %H:%M:%S",tz='UTC')  
  }
  diffs<-as.numeric(c(difftime(fulltimes,fulltimes[1],units="secs")))
  diffs<-diffs+(ms/1000)
  diffs
}  
CalculateBaseline=function(dfm){
  window.min=dfm$Parameters$Baseline.Window.Minutes
  newData<-dfm$RawData
  # the number of samples in those minutes
  window<-window.min*60*5
  if(window %% 2 ==0) 
    window=window+1
  
  for(i in 1:12) {
    cname <-paste("W",i,sep="")
    tmp<-runmed(newData[,cname],window)
    newData[,cname]<-newData[,cname]-tmp
  }
  
  dfm$BaselineData=newData
  ## Everything else must be recalculated
  dfm<-SetThreshold(dfm)
  dfm
}
SetThreshold = function(dfm,getStandard=TRUE) {
  ## First set the threshold...
  if(is.null(dfm$BaselineData)) {
    stop("DFM must have baseline.")
  }
  
  dfm<-Set.Fixed.Threshold(dfm)    
  
  ## Now update the licks and PI
  dfm<-Set.Feeding.Data(dfm)
  dfm<-Set.Tasting.Data(dfm)
  
  #Other measures
  dfm<-Set.Durations.And.Intervals(dfm)
  dfm<-Set.Tasting.Durations.And.Intervals(dfm)
  dfm
  
}
Set.Feeding.Data<-function(dfm){
  if(is.null(dfm$BaselineData))
    stop("Baseline must be calculated")
  
  newData<-dfm$BaselineData
  newData2<-dfm$BaselineData
  for(i in 1:12) {
    tmp<-Set.Feeding.Data.Well(dfm,i)
    cname <-paste("W",i,sep="")
    newData[,cname]<-tmp[,1]
    newData2[,cname]<-tmp[,2]
  }
  dfm$LickData<-newData
  dfm$EventData<-newData2
  dfm  
}
Set.Feeding.Data.Well<-function(dfm,well){
  ## Get all possible feeding Licks
  thresh<-Thresholds.Well(dfm,well)
  data<-BaselinedData.Well(dfm,well)
  
  Feeding.Licks.Min<-(data > thresh$FeedingMin)
  
  Feeding.Licks.Max<-(data > thresh$FeedingMax)
  
  ## Find continguous events above min threshold with at least one value above max threshold.
  ## The result of this function is also equivalent to the Events vector
  Events<-Get.Surviving.Events(Feeding.Licks.Min,Feeding.Licks.Max)
  
  ## Now remove events that are too short
  Events[Events<dfm$Parameters$Feeding.Minevents]<-0
  
  ## Now expand the licks to TRUE/FALSE entries
  FeedingLicks<-Expand.Events(Events)
  
  ## Now Bridge sporadic lick events into single events.
  tmp<-Link.Events(FeedingLicks,dfm$Parameters$Feeding.Event.Link.Gap)
  Events<-Get.Events(tmp)
  
  data.frame(FeedingLicks,Events)
}
Set.Tasting.Data<-function(dfm){
  if(is.null(dfm$BaselineData))
    stop("Baseline must be calculated")
  if(is.null(dfm$LickData))
    stop("Feeding Licks must be calculated")
  
  newData<-dfm$BaselineData
  newData2<-dfm$BaselineData
  
  for(i in 1:12) {
    tmp<-Set.Tasting.Data.Well(dfm,i)
    cname <-paste("W",i,sep="")
    newData[,cname]<-tmp[,1]
    newData2[,cname]<-tmp[,2]
  }
  dfm$TastingData<-newData
  dfm$TastingEventData<-newData2
  dfm
}
Set.Tasting.Data.Well<-function(dfm,well){
  ## Get Tasting Licks
  ## Note that Feeding Licks have to be calculated first because if the fly is 
  ## feeding, then tasting events have to be cancelled.
  thresh<-Thresholds.Well(dfm,well)
  data<-BaselinedData.Well(dfm,well)
  
  Licks<-(data > thresh$TastingMin & 
            data < thresh$TastingMax)
  
  FeedingLicks<-FeedingData.Well.Licks(dfm,well)
  
  ## Keep only taste licks that are not feeding licks
  Licks[FeedingLicks]<-FALSE
  
  Events<-Get.Events(Licks)
  
  ## Now remove events that are too short
  Events[Events<dfm$Parameters$Tasting.Minevents]<-0
  
  data.frame(Licks,Events)
}
Set.Fixed.Threshold<-function(dfm){
  tmp<-Set.Fixed.Threshold.Well(dfm,1)
  Thresholds=list(W1=tmp)
  
  for(i in 2:12){
    s<-paste("W",i,sep="")
    tmp<-Set.Fixed.Threshold.Well(dfm,i)
    Thresholds[[s]]<-tmp    
  }
  dfm$Thresholds<-Thresholds
  dfm  
}
Set.Fixed.Threshold.Well<-function(dfm,well){
  n<-SampleCount(dfm)
  ## Get well specific thresholds if the values are < 0 
  if(dfm$Parameters$Feeding.Threshold<0){
    ## Find maximum reading
    tmp<-max(BaselinedData.Well(dfm,well))
    tmpA <- round(tmp*abs(dfm$Parameters$Feeding.Threshold),0)
    tmpB <- round(tmp*abs(dfm$Parameters$Feeding.Minimum),0)
    tmpC <- round(tmp*abs(dfm$Parameters$Tasting.Minimum),0)
    tmpD <-round(tmp*abs(dfm$Parameters$Tasting.Maximum),0)
  }
  else {
    tmpA<-dfm$Parameters$Feeding.Threshold 
    tmpB<-dfm$Parameters$Feeding.Minimum
    tmpC<-dfm$Parameters$Tasting.Minimum
    tmpD<-dfm$Parameters$Tasting.Maximum
  }
  
  feeding.max.thresh<-rep(tmpA,n)
  feeding.min.thresh<-rep(tmpB,n)
  tasting.min.thresh<-rep(tmpC,n)
  tasting.max.thresh<-rep(tmpD,n)
  
  
  r.tmp<-data.frame(feeding.max.thresh,feeding.min.thresh,tasting.max.thresh,tasting.min.thresh)       
  
  names(r.tmp)<-c("FeedingMax","FeedingMin","TastingMax","TastingMin")                            
  r.tmp  
}
Set.Durations.And.Intervals<-function(dfm){
  tmp<-Set.Durations.And.Intervals.Well(dfm,1)
  Durations = list(W1=tmp$Durations)
  Intervals = list(W1=tmp$Intervals)
  
  for(i in 2:12){
    s<-paste("W",i,sep="")
    tmp<-Set.Durations.And.Intervals.Well(dfm,i)
    Durations[[s]]<-tmp$Durations  
    Intervals[[s]]<-tmp$Intervals
  }
  
  dfm$Durations<-Durations
  dfm$Intervals<-Intervals
  dfm
}
Set.Durations.And.Intervals.Well<-function(dfm,well){
  data<-BaselineData.Well(dfm,well)
  events<-FeedingData.Well.Events(dfm,well)
  ## Now we need to update the event durations
  ## Indices will be used for summary duration characteristics
  indices<-1:length(events)
  
  indices<-indices[events>0]
  boutDurs<-events[events>0]  
  
  Durations<-0
  
  if(length(boutDurs)>0) {
    max.inten<-rep(0,length(indices))
    min.inten<-rep(0,length(indices))
    sum.inten<-rep(0,length(indices))
    avg.inten<-rep(0,length(indices))    
    var.inten<-rep(0,length(indices))   
    for(i in 1:length(indices)){
      dataindex<-indices[i]
      eventlength<-boutDurs[i]
      tmp2<-data[dataindex:(dataindex+(eventlength-1))]
      max.inten[i]<-max(tmp2)
      min.inten[i]<-min(tmp2)
      sum.inten[i]<-sum(tmp2)
      avg.inten[i]<-mean(tmp2)  
      var.inten[i]<-var(tmp2)
    }
    
    BoutData<-data.frame(min.inten,max.inten,sum.inten,avg.inten,var.inten)
    names(BoutData)<-c("MinIntensity","MaxIntensity","SumIntensity","MeanIntensity","VarIntensity")
    
    tmp<-BaselineData(dfm)
    tmp<-tmp[indices,]
    Minutes<-tmp$Minutes
    Events<-boutDurs
    Duration<-Events/dfm$Parameters$Samples.Per.Sec
    AvgInten<-BoutData$MeanIntensity
    MaxInten<-BoutData$MaxIntensity
    MinInten<-BoutData$MinIntensity
    SumInten<-BoutData$SumIntensity
    VarInten<-BoutData$VarIntensity
    Durations<-data.frame(Minutes,Events,Duration,SumInten,AvgInten,MinInten,MaxInten,VarInten)
    names(Durations)<-c("Minutes","Licks","Duration","TotalIntensity","AvgIntensity","MinIntensity","MaxIntensity","VarIntensity")    
  }
  
  result<-list(Durations=Durations)
  
  ## Now intervals
  
  ## Collapse feeding data to time BETWEEN events.
  ##boutInt<-Get.Intervals(FeedingData.Well.Licks(dfm,well))  
  tmp<-FeedingData.Well.Events(dfm,well)
  tmp<-Expand.Events(tmp)
  boutInt<-Get.Intervals(tmp)  
  
  
  indices<-1:length(boutInt)
  indices<-indices[boutInt>0]
  boutInt<-boutInt[boutInt>0]
  
  
  spm<-dfm$Parameters$Samples.Per.Sec
  intA<-boutInt/spm
  
  Ints<-0
  
  if(length(intA)>0) {
    tmp<-BaselineData(dfm)
    tmp<-tmp[indices,]
    Minutes<-tmp$Minutes
    Sample<-tmp$Sample
    IntervalSec<-intA
    Ints<-data.frame(Minutes,Sample,IntervalSec)
  }
  
  result<-list(Durations=Durations,Intervals=Ints)
  result
}
Set.Tasting.Durations.And.Intervals<-function(dfm){
  tmp<-Set.Tasting.Durations.And.Intervals.Well(dfm,1)
  Durations = list(W1=tmp$Durations)
  Intervals = list(W1=tmp$Intervals)
  
  for(i in 2:12){
    s<-paste("W",i,sep="")
    tmp<-Set.Tasting.Durations.And.Intervals.Well(dfm,i)
    Durations[[s]]<-tmp$Durations  
    Intervals[[s]]<-tmp$Intervals
  }
  
  dfm$TastingDurations<-Durations
  dfm$TastingIntervals<-Intervals
  dfm
}
Set.Tasting.Durations.And.Intervals.Well<-function(dfm,well){
  data<-BaselineData.Well(dfm,well)
  events<-TastingData.Well.Events(dfm,well)
  ## Now we need to update the event durations
  ## Indices will be used for summary duration characteristics
  indices<-1:length(events)
  
  indices<-indices[events>0]
  boutDurs<-events[events>0]  
  
  Durations<-0
  
  if(length(boutDurs)>0) {
    max.inten<-rep(0,length(indices))
    min.inten<-rep(0,length(indices))
    sum.inten<-rep(0,length(indices))
    avg.inten<-rep(0,length(indices))    
    var.inten<-rep(0,length(indices))   
    for(i in 1:length(indices)){
      dataindex<-indices[i]
      eventlength<-boutDurs[i]
      tmp2<-data[dataindex:(dataindex+(eventlength-1))]
      max.inten[i]<-max(tmp2)
      min.inten[i]<-min(tmp2)
      sum.inten[i]<-sum(tmp2)
      avg.inten[i]<-mean(tmp2)  
      var.inten[i]<-var(tmp2)
    }
    
    BoutData<-data.frame(min.inten,max.inten,sum.inten,avg.inten,var.inten)
    names(BoutData)<-c("MinIntensity","MaxIntensity","SumIntensity","MeanIntensity","VarIntensity")
    
    tmp<-BaselineData(dfm)
    tmp<-tmp[indices,]
    Minutes<-tmp$Minutes
    Events<-boutDurs
    Duration<-Events/dfm$Parameters$Samples.Per.Sec
    AvgInten<-BoutData$MeanIntensity
    MaxInten<-BoutData$MaxIntensity
    MinInten<-BoutData$MinIntensity
    SumInten<-BoutData$SumIntensity
    VarInten<-BoutData$VarIntensity
    Durations<-data.frame(Minutes,Events,Duration,SumInten,AvgInten,MinInten,MaxInten,VarInten)
    names(Durations)<-c("Minutes","Licks","Duration","TotalIntensity","AvgIntensity","MinIntensity","MaxIntensity","VarIntensity")    
  }
  
  result<-list(Durations=Durations)
  
  ## Now intervals
  
  ## Collapse feeding data to time BETWEEN events.
  boutInt<-Get.Intervals(FeedingData.Well.Licks(dfm,well))  
  
  indices<-1:length(boutInt)
  indices<-indices[boutInt>0]
  boutInt<-boutInt[boutInt>0]
  
  
  spm<-dfm$Parameters$Samples.Per.Sec
  intA<-boutInt/spm
  
  Ints<-0
  
  if(length(intA)>0) {
    tmp<-BaselineData(dfm)
    tmp<-tmp[indices,]
    Minutes<-tmp$Minutes
    Sample<-tmp$Sample
    IntervalSec<-intA
    Ints<-data.frame(Minutes,Sample,IntervalSec)
  }
  
  result<-list(Durations=Durations,Intervals=Ints)
  result
}
Thresholds.Well<-function(dfm,well,range=c(0,0)){
  cname=paste("W",well,sep="")
  tmp<-dfm$Thresholds[[cname]]
  if(sum(range)!=0) {
    tmp<- tmp[(dfm$BaselineData$Minutes>range[1]) & (dfm$BaselineData$Minutes<=range[2]),]
  }    
  tmp  
}

##### Data Access functions #####
BaselinedData.Well<-function(dfm,well,range=c(0,0)) {  
  cname=paste("W",well,sep="")
  tmp<-dfm$BaselineData[,cname]  
  if(sum(range)!=0) {
    tmp<- tmp[(dfm$BaselineData$Minutes>range[1]) & (dfm$BaselineData$Minutes<=range[2])]
  }    
  tmp  
}
BaselinedData<-function(dfm,range=c(0,0)) {   
  tmp<-dfm$BaselineData
  if(sum(range)!=0) {
    tmp<- tmp[(dfm$BaselineData$Minutes>range[1]) & (dfm$BaselineData$Minutes<=range[2]),]
  }    
  tmp  
}
SampleCount<-function(dfm,range=c(0,0)){
  nrow(BaselinedData(dfm,range))  
}
FeedingData.Well.Licks<-function(dfm,well,range=c(0,0)){
  cname=paste("W",well,sep="")
  tmp<-FeedingData.Licks(dfm,range)
  tmp[,cname]  
}
## Remember that this function returns a vector with 
## duration of event information as well.
## Need to set these to 1 to get number of events.
FeedingData.Well.Events<-function(dfm,well,range=c(0,0)){
  cname=paste("W",well,sep="")
  tmp<-FeedingData.Events(dfm,range)
  tmp[,cname]  
}
TastingData.Well<-function(dfm,well,range=c(0,0)){
  cname=paste("W",well,sep="")
  tmp<-dfm$TastingData[,cname]  
  if(sum(range)!=0) {
    tmp<- tmp[(tmp$Minutes>range[1]) & (tmp$Minutes<=range[2])]
  }   
  tmp    
}
TastingData.Well.Events<-function(dfm,well,range=c(0,0)){
  cname=paste("W",well,sep="")
  tmp<-TastingData.Events(dfm,range)
  tmp[,cname]  
}
TastingData.Events<-function(dfm,range=c(0,0)){
  data<-dfm$TastingEventData  
  if(sum(range)!=0) {
    data<- data[(data$Minutes>range[1] & data$Minutes<=range[2]),]
  }    
  data
}
FeedingData.Licks<-function(dfm,range=c(0,0)){
  data<-dfm$LickData  
  if(sum(range)!=0) {
    data<- data[(data$Minutes>range[1] & data$Minutes<=range[2]),]
  }    
  data
}
## Remember that this function returns a vector with 
## duration of event information as well.
## Need to set these to 1 to get number of events.
FeedingData.Events<-function(dfm,range=c(0,0)){
  data<-dfm$EventData  
  if(sum(range)!=0) {
    data<- data[(data$Minutes>range[1] & data$Minutes<=range[2]),]
  }    
  data
}
TastingData<-function(dfm,range=c(0,0)){
  data<-dfm$TastingData  
  if(sum(range)!=0) {
    data<- data[(data$Minutes>range[1] & data$Minutes<=range[2]),]
  }    
  data
}
Feeding.TotalLicks<-function(dfm,range=c(0,0)){
  result<-rep(-1,12)
  data<-FeedingData.Licks(dfm,range)
  for(i in 1:12) {
    cname=paste("W",i,sep="")
    tmp<-data[,cname]
    result[i]<-sum(tmp)
  }
  names(result)<-paste("W",1:12,sep="")
  result
}
Feeding.TotalLicks.Well<-function(dfm,well,range=c(0,0)){
  tmp<-Feeding.TotalLicks(dfm,range)
  tmp[well]
}
Feeding.TotalEvents<-function(dfm,range=c(0,0)){
  result<-rep(-1,12)
  data<-FeedingData.Events(dfm,range)
  for(i in 1:12) {
    cname=paste("W",i,sep="")
    tmp<-data[,cname]
    result[i]<-sum(tmp>0)
  }
  names(result)<-paste("W",1:12,sep="")
  result  
}
Feeding.TotalEvents.Well<-function(dfm,well,range=c(0,0)){
  tmp<-Feeding.TotalEvents(dfm,range)
  tmp[well]
}
Tasting.TotalLicks<-function(dfm,range=c(0,0)){
  result<-rep(-1,12)
  data<-TastingData(dfm,range)
  for(i in 1:12) {
    cname=paste("W",i,sep="")
    tmp<-data[,cname]
    result[i]<-sum(tmp)
  }
  names(result)<-paste("W",1:12,sep="")
  result
}
Tasting.TotalLicks.Well<-function(dfm,well,range=c(0,0)){
  tmp<-Tasting.TotalLicks(dfm,range)
  tmp[well]
}
BaselineData.Well=function(dfm,well,range=c(0,0)) {
  cname=paste("W",well,sep="")
  tmp<-BaselineData(dfm,range)
  tmp[,cname]
}
BaselinedData.Range.Well<-function(dfm,well,range=c(0,0)){
  tmp<-BaselinedData.Well(dfm,well,range)
  x1<-min(tmp)
  x2<-max(tmp)
  c(x1,x2)
}
Minutes<-function(dfm,range=c(0,0)) {
  data<-dfm$BaselineData$Minutes
  if(sum(range)!=0) {
    data<- data[(data>range[1] & data<range[2])]
  }    
  data
}
Feeding.Durations.Well<-function(dfm,well){
  cname=paste("W",well,sep="")
  adurs<-dfm$Durations[[cname]]
  adurs
}
Feeding.Intervals.Well<-function(dfm,well){
  cname=paste("W",well,sep="")
  adurs<-dfm$Intervals[[cname]]
  adurs
}
LastSampleData.Well<-function(dfm,well){
  tmp<-BaselinedData.Well(dfm,well)
  tmp[length(tmp)]
}
FirstSampleData.Well<-function(dfm,well){
  tmp<-BaselinedData.Well(dfm,well)  
  tmp[1]
}
## These functions will output the intervals
GetIntervalData.Well<-function(dfm,well, range=c(0,0)){
  nameA<-paste("W",well,sep="")
  parameter.vector<-GetDFMParameterVector(dfm)
  pnames<-Get.Parameter.Names(dfm$Parameters)
  
  theData<-dfm$Intervals[[nameA]]
  if(length(theData)==1 && theData[1]==0){
    theData<-data.frame(matrix(rep(NA,8),nrow=1))
  }
  else {
    if(sum(range)!=0) {
      theData<- theData[(theData$Minutes>range[1] & theData$Minutes<=range[2]),]
    }
  }
  
  Well<-rep(well,nrow(theData))
  chamber<-rep(GetChamberFromWell(dfm,well),nrow(theData))
  TCWell<-rep(GetTCWellFromWell(dfm,well),nrow(theData))
  DFM<-rep(dfm$ID,nrow(theData))
  
  tmpA<-data.frame(DFM,chamber,TCWell,Well,theData)
  tmp2<-matrix(rep(parameter.vector,nrow(tmpA)),ncol=length(parameter.vector),byrow=TRUE)
  tmp3<-data.frame(tmpA,tmp2)
  names(tmp3)<-c("DFM","Chamber","TCWell","Well","Minutes","Sample","IntervalSec",pnames)
  
  if(dfm$Parameters$Chamber.Size==1)
    tmp3<-tmp3[,!names(tmp3) %in% "TCWell"]
  
  tmp3
}
GetDurationData.Well<-function(dfm,well, range=c(0,0)){
  nameA<-paste("W",well,sep="")
  parameter.vector<-GetDFMParameterVector(dfm)
  pnames<-Get.Parameter.Names(dfm$Parameters)
  
  theData<-dfm$Durations[[nameA]]
  if(length(theData)==1 && theData[1]==0){
    theData<-data.frame(matrix(rep(NA,8),nrow=1))
  }
  else {
    if(sum(range)!=0) {
      theData<- theData[(theData$Minutes>range[1] & theData$Minutes<=range[2]),]
    }
    if(nrow(theData)==0){
      theData<-data.frame(matrix(rep(NA,8),nrow=1))
    }
  }
  
  
  
  Well<-rep(well,nrow(theData))
  chamber<-rep(GetChamberFromWell(dfm,well),nrow(theData))
  TCWell<-rep(GetTCWellFromWell(dfm,well),nrow(theData))
  DFM<-rep(dfm$ID,nrow(theData))
  
  tmpA<-data.frame(DFM,chamber,TCWell,Well,theData)
  tmp2<-matrix(rep(parameter.vector,nrow(tmpA)),ncol=length(parameter.vector),byrow=TRUE)
  tmp3<-data.frame(tmpA,tmp2)
  names(tmp3)<-c("DFM","Chamber","TCWell","Well","Minutes","Licks","Duration","TotalIntensity","AvgIntensity","MinIntensity","MaxIntensity","VarIntensity",pnames)
  
  if(dfm$Parameters$Chamber.Size==1)
    tmp3<-tmp3[,!names(tmp3) %in% "TCWell"]
  
  tmp3
}
GetTCWellFromWell<-function(dfm,well){
  if(dfm$Parameters$Chamber.Size==1){
    well<-NA
  }
  else {
    if(well==1 || well==3 || well==5 || well==7 || well==9 || well==11){
      if(dfm$Parameters$PI.Multiplier==1) {
        well<-"WellA"
      }  
      else {
        well<-"WellB"
      }
    }
    
    else {
      if(dfm$Parameters$PI.Multiplier==1) {
        well<-"WellB"
      }  
      else {
        well<-"WellA"
      } 
    }
  }
  well
}
GetChamberFromWell<-function(dfm,well){
  if(dfm$Parameters$Chamber.Size==1){
    chamber<-well
  }
  else {
    if(well==1 || well==2)
      chamber<-1
    else if(well==3 || well==4)
      chamber<-2
    else if(well==5 || well==6)
      chamber<-3
    else if(well==7 || well==8)
      chamber<-4
    else if(well==9 || well==10)
      chamber<-5
    else if(well==11 || well==12)
      chamber<-6
  }
  chamber
}
UpdateHiddenDFMObject<-function(dfm){
  st<-paste("DFM",dfm$ID,sep="")
  assign(st,dfm,pos=1)
}
GetDFMParameterVector<-function(dfm){
  GetParameterVector(dfm$Parameters)  
}


##### Analysis helper functions #####
Feeding.IntervalSummary.Well<-function(dfm,well,range=c(0,0)){
  cname=paste("W",well,sep="")
  adurs<-dfm$Intervals[[cname]]
  if(sum(range)!=0){
    if(!is.data.frame(adurs)){
      a<-0
      aa<-0
    }
    else {
      adurs<-adurs[(adurs$Minutes>range[1]) & (adurs$Minutes<=range[2]),]  
      if(nrow(adurs)==0){
        a<-0
        aa<-0
      }
      else {    
        a<-mean(adurs$IntervalSec)
        aa<-median(adurs$IntervalSec)  
      }
    }  
  }
  else {
    if(!is.data.frame(adurs)){
      a<-0
      aa<-0
    } else {
      a<-mean(adurs$IntervalSec)  
      aa<-median(adurs$IntervalSec)
    }
  }
  
  if(is.na(a)||is.nan(a)) a<-0
  if(is.na(aa)||is.nan(aa)) aa<-0
  tmp<-data.frame(a,aa)
  names(tmp)<-c("MeanTimeBtw","MedTimeBtw")
  tmp
}
Feeding.DurationSummary.Well<-function(dfm,well,range=c(0,0)){
  cname=paste("W",well,sep="")
  adurs<-dfm$Durations[[cname]]
  if(sum(range)!=0){
    if(!is.data.frame(adurs)){
      a<-NA
      aa<-NA
    }
    else {
      adurs<-adurs[(adurs$Minutes>range[1]) & (adurs$Minutes<=range[2]),]  
      if(nrow(adurs)==0){
        a<-NA
        aa<-NA
      }
      else {    
        a<-mean(adurs$Duration)
        aa<-median(adurs$Duration)  
      }
    }
  }
  else {
    if(!is.data.frame(adurs)){
      a<-NA
      aa<-NA
    } else {
      a<-mean(adurs$Duration)  
      aa<-median(adurs$Duration)
    }
  }
  
  if(is.na(a)||is.nan(a)) a<-NA
  if(is.na(aa)||is.nan(aa)) aa<-NA
  tmp<-data.frame(a,aa)
  
  names(tmp)<-c("MeanDur","MedianDur")
  tmp
}
Feeding.IntensitySummary.Well<-function(dfm,well,range=c(0,0)){
  d<-BaselineData.Well(dfm,well,range)
  l<-Expand.Events(FeedingData.Well.Events(dfm,well,range))
  
  da<-d[l]
  
  if(length(da)==0){
    a<-NA
    aa<-NA
    aaa<-NA
    aaaa<-NA
  }
  else {
    a<-mean(da)  
    aa<-median(da)
    aaa<-min(da)
    aaaa<-max(da)
  }
  
  tmp<-data.frame(a,aa,aaa,aaaa)
  names(tmp)<-c("MeanInt","MedianInt","MinInt","MaxInt")
  tmp
}
Feeding.Summary.OneWell<-function(dfm,range=c(0,0),TransformLicks=TRUE){
  if(dfm$Parameters$Chamber.Size!=1)
    stop("This function is for single chambers only")
  lights.sec<-(apply(GetLightsInfo(dfm,range)[,2:13],2,sum))/dfm$Parameters$Samples.Per.Second
  for(i in 1:12 ){
    interval<-Feeding.IntervalSummary.Well(dfm,i,range)
    intensity<-Feeding.IntensitySummary.Well(dfm,i,range)
    dur<-Feeding.DurationSummary.Well(dfm,i,range)  
    FLicks<-Feeding.TotalLicks.Well(dfm,i,range)
    FEvents<-Feeding.TotalEvents.Well(dfm,i,range)    
    if(i==1)
      result<-data.frame(matrix(c(dfm$ID,i,FLicks,FEvents,unlist(dur),unlist(interval),unlist(intensity),lights.sec[i],range[1],range[2]),nrow=1))  
    else {
      tmp<-data.frame(matrix(c(dfm$ID,i,FLicks,FEvents,unlist(dur),unlist(interval),unlist(intensity),lights.sec[i],range[1],range[2]),nrow=1))  
      result<-rbind(result,tmp)
    }      
  }
  
  names(result)<-c("DFM","Chamber","Licks","Events","MeanDuration","MedDuration",
                   "MeanTimeBtw","MedTimeBtw","MeanInt","MedianInt","MinInt","MaxInt","OptoOn_sec","StartMin","EndMin")
  if(TransformLicks==TRUE)
    result$Licks<-result$Licks^0.25
  result    
  
}
Feeding.Summary.TwoWell<-function(dfm,range=c(0,0),TransformLicks=TRUE){
  if(dfm$Parameters$Chamber.Size!=2)
    stop("This function is for two-chamber DFM only")
  
  lights.sec<-(apply(GetLightsInfo(dfm,range)[,2:13],2,sum))/dfm$Parameters$Samples.Per.Second
  
  for(i in 1:nrow(dfm$Parameters$Chamber.Sets)) {
    if(dfm$Parameters$PI.Multiplier==1){
      wellA<-dfm$Parameters$Chamber.Sets[i,1]
      wellB<-dfm$Parameters$Chamber.Sets[i,2] 
    }
    else {
      wellB<-dfm$Parameters$Chamber.Sets[i,1]
      wellA<-dfm$Parameters$Chamber.Sets[i,2] 
    }
    
    interval.a<-Feeding.IntervalSummary.Well(dfm,wellA,range)
    intensity.a<-Feeding.IntensitySummary.Well(dfm,wellA,range)
    dur.a<-Feeding.DurationSummary.Well(dfm,wellA,range)  
    FLicks.a<-Feeding.TotalLicks.Well(dfm,wellA,range)
    FEvents.a<-Feeding.TotalEvents.Well(dfm,wellA,range)    
    
    interval.b<-Feeding.IntervalSummary.Well(dfm,wellB,range)
    intensity.b<-Feeding.IntensitySummary.Well(dfm,wellB,range)
    dur.b<-Feeding.DurationSummary.Well(dfm,wellB,range)  
    FLicks.b<-Feeding.TotalLicks.Well(dfm,wellB,range)
    FEvents.b<-Feeding.TotalEvents.Well(dfm,wellB,range) 
    
    lights.a<-lights.sec[wellA]
    lights.b<-lights.sec[wellB]
    
    FPIs<-c((FLicks.a-FLicks.b)/(FLicks.a+FLicks.b),(FEvents.a-FEvents.b)/(FEvents.a+FEvents.b))
    
    
    
    if(i==1){
      result<-data.frame(matrix(c(dfm$ID,i,FPIs,FLicks.a,FLicks.b,FEvents.a,FEvents.b,unlist(dur.a),unlist(dur.b),unlist(interval.a),
                                  unlist(interval.b),unlist(intensity.a),unlist(intensity.b),lights.a,lights.b,range[1],range[2]),nrow=1))    
      
    }
    else {
      tmp<-data.frame(matrix(c(dfm$ID,i,FPIs,FLicks.a,FLicks.b,FEvents.a,FEvents.b,unlist(dur.a),unlist(dur.b),unlist(interval.a),
                               unlist(interval.b),unlist(intensity.a),unlist(intensity.b),lights.a,lights.b,range[1],range[2]),nrow=1))
      result<-rbind(result,tmp)      
    }
  }
  names(result)<-c("DFM","Chamber","PI","EventPI","LicksA","LicksB","EventsA","EventsB","MeanDurationA","MedDurationA",
                   "MeanDurationB","MedDurationB","MeanTimeBtwA","MedTimeBtwA",
                   "MeanTimeBtwB","MedTimeBtwB","MeanIntA","MedianIntA","MinIntA","MaxIntA",
                   "MeanIntB","MedianIntB","MinIntB","MaxIntB","OptoOn_sec_A","OptoOn_sec_B","StartMin","EndMin")
  if(TransformLicks==TRUE){
    result$LicksA<-result$LicksA^0.25
    result$LicksB<-result$LicksB^0.25
  }
  result    
  
}
AggregateTreatmentsBinnedData<-function(results){
  trt.summary1<-aggregate(results,by=list(results$Interval,results$Treatment),mean,na.rm=TRUE) 
  trt.summary2<-aggregate(results,by=list(results$Interval,results$Treatment),mySEM)
  trt.summary2<-trt.summary2[,-grep("Treatment|DFM|Chamber|Interval",colnames(trt.summary2))]
  trt.summary1<-trt.summary1[,-grep("Treatment|DFM|Chamber|Interval",colnames(trt.summary1))]
  
  ## This just indicates a two well chamber
  if("LicksA" %in% names(results)){
    tmp<-names(trt.summary1)[4:25]
    tmps<-paste(tmp,"SEM",sep="")
    tmp<-names(trt.summary1)
    tmp<-c(tmp[1:25],tmps,tmp[26:ncol(trt.summary1)])
    tmp[1]<-"Interval"
    tmp[2]<-"Treatment"
    trt.summary<-data.frame(trt.summary1[,1:25],trt.summary2[,4:25],trt.summary1[,26:ncol(trt.summary1)])
    names(trt.summary)<-tmp
  }
  else {
    trt.summary<-data.frame(trt.summary1[,1:13],trt.summary2[,4:13],trt.summary1[,14:ncol(trt.summary1)])
    names(trt.summary1)[names(trt.summary1) == "Group.1"] <- "Interval"
    names(trt.summary1)[names(trt.summary1) == "Group.2"] <- "Treatment"
    tmp<-names(trt.summary1)[4:13]
    tmp<-paste(tmp,"SEM",sep="")
    names(trt.summary)<-c(names(trt.summary1)[1:13],tmp,names(trt.summary1)[14:ncol(trt.summary1)])
  }
  trt.summary
}
AggregateTreatments<-function(results){
  trt.summary1<-aggregate(results,by=list(results$Treatment),mean, na.rm=TRUE) 
  trt.summary2<-aggregate(results,by=list(results$Treatment),mySEM)
  trt.summary1<-trt.summary1[,-grep("Treatment|DFM|Chamber",colnames(trt.summary1))]
  trt.summary2<-trt.summary2[,-grep("Treatment|DFM|Chamber",colnames(trt.summary2))]
  
  if("LicksA" %in% names(results)){
    names(trt.summary1)[names(trt.summary1) == "Group.1"] <- "Treatment"
    tmp<-names(trt.summary1)[2:23]
    tmp<-paste(tmp,"SEM",sep="")
    trt.summary<-data.frame(trt.summary1[,1:23],trt.summary2[,2:23],trt.summary1[,24:ncol(trt.summary1)])
    names(trt.summary)<-c(names(trt.summary1)[1:23],tmp,names(trt.summary1)[24:ncol(trt.summary1)])
  }
  else {
    names(trt.summary1)[names(trt.summary1) == "Group.1"] <- "Treatment"
    tmp<-names(trt.summary1)[2:11]
    tmp<-paste(tmp,"SEM",sep="")
    trt.summary<-data.frame(trt.summary1[,1:11],trt.summary2[,2:11],trt.summary1[,12:ncol(trt.summary1)])
    names(trt.summary)<-c(names(trt.summary1)[1:11],tmp,names(trt.summary1)[12:ncol(trt.summary1)])
  }
  
  
  #tmp<-c("Treatment","MeanLicks","MeanEvents","MeanMDuration","MeanMedDuration","MeanMTimeBtw","MeanMedTimeBtw","MeanMInt","MeanMedInt",
  #       "SEMLicks","SEMEvents","SEMMDuration","SEMMedDuration","SEMMTimeBtw","SEMMedTimeBtw","SEMMInt","SEMMedInt",names(trt.summary)[18:ncol(trt.summary)])
  
  #names(trt.summary)<-tmp
  trt.summary
}


##### Data output helper functions. These can not be generalized. #####
OutputBaselinedData.Monitors<-function(monitors,parameters,range=c(0,0),filename="Baselined"){
  individ.params<-FALSE
  ## Check to determine whether parameters is a signle parameter object
  ## or a list of them.  If it is a single one, then we use the same one for all
  ## if it is a list, then we use a different one for each.
  if(is.list(parameters[[1]])==TRUE){
    if(length(parameters)!=length(monitors))
      stop("If individuals parameter objects are specified, there must be one for each DFM.")
    individ.params<-TRUE
  }
  
  for(j in 1:length(monitors)){
    print(paste("Outputting Baselined Data for DFM ",monitors[j],".",sep=""))
    flush.console()
    monitor<-monitors[j]
    if(individ.params==TRUE)
      p<-parameters[[j]]
    else
      p<-parameters
    dfm<-DFMClass(monitor,p)
    OutputBaselinedData.DFM(dfm,range,filename)
  }
}
OutputIntervalData.Monitors<-function(monitors,parameters,expDesign=NA,range=c(0,0),filename="IntervalData"){
  individ.params<-FALSE
  ## Check to determine whether parameters is a signle parameter object
  ## or a list of them.  If it is a single one, then we use the same one for all
  ## if it is a list, then we use a different one for each.
  if(is.list(parameters[[1]])==TRUE){
    if(length(parameters)!=length(monitors))
      stop("If individuals parameter objects are specified, there must be one for each DFM.")
    individ.params<-TRUE
  }
  for(j in 1:length(monitors)){
    ##print(paste("Outputting Interval Data for DFM ",monitors[j],".",sep=""))
    ##flush.console()
    monitor<-monitors[j]
    if(individ.params==TRUE)
      p<-parameters[[j]]
    else
      p<-parameters
    dfm<-DFMClass(monitor,p)
    tmp2<-GetIntervalData.DFM(dfm,range)
    if(is.data.frame(expDesign))
      tmp2<-AppendTreatmentonResultsFrame(tmp2,expDesign)
    if(j==1){
      result<-tmp2
    }
    else {
      result<-rbind(result,tmp2)
    }
  }
  filename<-paste(filename,".csv",sep="") 
  write.csv(result,file=filename,row.names=FALSE)  
}
OutputDurationData.Monitors<-function(monitors,parameters,expDesign=NA,range=c(0,0),filename="DurationsData"){
  individ.params<-FALSE
  ## Check to determine whether parameters is a signle parameter object
  ## or a list of them.  If it is a single one, then we use the same one for all
  ## if it is a list, then we use a different one for each.
  if(is.list(parameters[[1]])==TRUE){
    if(length(parameters)!=length(monitors))
      stop("If individuals parameter objects are specified, there must be one for each DFM.")
    individ.params<-TRUE
  }
  for(j in 1:length(monitors)){
    ##print(paste("Outputting Interval Data for DFM ",monitors[j],".",sep=""))
    ##flush.console()
    monitor<-monitors[j]
    if(individ.params==TRUE)
      p<-parameters[[j]]
    else
      p<-parameters
    dfm<-DFMClass(monitor,p)
    tmp2<-GetDurationData.DFM(dfm,range)
    if(is.data.frame(expDesign))
      tmp2<-AppendTreatmentonResultsFrame(tmp2,expDesign)
    if(j==1){
      result<-tmp2
    }
    else {
      result<-rbind(result,tmp2)
    }
  }
  filename<-paste(filename,".csv",sep="") 
  write.csv(result,file=filename,row.names=FALSE)  
}
## This fucntion will output for each well in each chamber for each monitor
## the total amount of time spend drinking over the perscribed range.
OutputTotalFeeding.Monitors<-function(monitors,parameters,expDesign=NA,range=c(0,0),filename="TotalFeedingTime"){
  individ.params<-FALSE
  ## Check to determine whether parameters is a signle parameter object
  ## or a list of them.  If it is a single one, then we use the same one for all
  ## if it is a list, then we use a different one for each.
  if(is.list(parameters[[1]])==TRUE){
    if(length(parameters)!=length(monitors))
      stop("If individuals parameter objects are specified, there must be one for each DFM.")
    individ.params<-TRUE
  }
  
  for(j in 1:length(monitors)){
    monitor<-monitors[j]
    x<-1:12
    if(individ.params==TRUE)
      p<-parameters[[j]]
    else
      p<-parameters
    dfm<-DFMClass(monitor,p)
    parameter.vector<-GetParameterVector(p)
    pnames<-Get.Parameter.Names(p)
    tmp<-Feeding.Summary.DFM(dfm,range)
    if(p$Chamber.Size==1){
      atotal<-tmp$Events*tmp$MeanDuration
      d<-tmp$DFM
      c<-tmp$Chamber
      tmp2<-matrix(rep(parameter.vector,length(d)),ncol=length(parameter.vector),byrow=TRUE)
      tmp3<-data.frame(d,c,atotal,tmp2)
      names(tmp3)<-c("DFM","Chamber","TotalSec",pnames)
      if(j==1){
        result<-tmp3
      }
      else {
        result<-rbind(result,tmp3)  
      }
    }
    else if(p$Chamber.Size==2){
      atotal<-tmp$EventsA*tmp$MeanDurationA
      btotal<-tmp$EventsB*tmp$MeanDurationB
      d<-tmp$DFM
      c<-tmp$Chamber
      tmp2<-matrix(rep(parameter.vector,length(d)),ncol=length(parameter.vector),byrow=TRUE)
      tmp3<-data.frame(d,c,atotal,btotal,tmp2)
      names(tmp3)<-c("DFM","Chamber","ATotalSec","BTotalSec",pnames)
      if(j==1){
        result<-tmp3
      }
      else {
        result<-rbind(result,tmp3)  
      }      
    }
    else 
      stop("Feeding Summary not implemented for this DFM type.")    
  }
  
  if(is.data.frame(expDesign))
    result<-AppendTreatmentonResultsFrame(result,expDesign)
  
  filename<-paste(filename,".csv",sep="") 
  write.csv(result,file=filename,row.names=FALSE)  
}


##### Plotting helper functions. #####
BinnedPlot.OneWell.Trt<-function(binnedDataResult,Type="Licks",SaveToFile=FALSE){
  tmp<-binnedDataResult$Stats
  
  pd <- position_dodge(5) # move them .05 to the left and right
  
  if(Type=="Licks") {
    ylabel<-"Licks"
    gp<-ggplot(tmp,aes(x=Minutes,y=Licks,color=Treatment,group=Treatment)) + 
      geom_errorbar(aes(ymin=Licks-LicksSEM, ymax=Licks+LicksSEM,color=Treatment), width=.1, position=pd) +
      geom_line(position=pd,size=1) +
      geom_point(position=pd, size=4, shape=21, fill="white") +xlab("Minutes") + ylab(ylabel)
    filename<-paste("BinnedLicksPlots.pdf",sep="")
    analysis<-data.frame(binnedDataResult$Results$Interval,binnedDataResult$Results$Treatment,binnedDataResult$Results$Licks)
    names(analysis)<-c("Interval","Treatment","Y")
  }
  else if(Type=="Events") {
    ylabel<-"Events"
    gp<-ggplot(tmp,aes(x=Minutes,y=Events,color=Treatment,group=Treatment)) + 
      geom_errorbar(aes(ymin=Events-EventsSEM, ymax=Events+EventsSEM,color=Treatment), width=.1, position=pd) +
      geom_line(position=pd,size=1) +
      geom_point(position=pd, size=4, shape=21, fill="white") +xlab("Minutes") + ylab("Events")
    filename<-paste("BinnedEventsPlots.pdf",sep="")
    analysis<-data.frame(binnedDataResult$Results$Interval,binnedDataResult$Results$Treatment,binnedDataResult$Results$Events)
    names(analysis)<-c("Interval","Treatment","Y")
  }
  else if(Type=="Durations") {
    gp<-ggplot(tmp,aes(x=Minutes,y=MeanDuration,color=Treatment,group=Treatment)) + 
      geom_errorbar(aes(ymin=MeanDuration-MeanDurationSEM, ymax=MeanDuration+MeanDurationSEM,color=Treatment), width=.1, position=pd) +
      geom_line(position=pd,size=1) +
      geom_point(position=pd, size=4, shape=21, fill="white") +xlab("Minutes") + ylab("Mean Event Duration")
    filename<-paste("BinnedDurationPlots.pdf",sep="")
    analysis<-data.frame(binnedDataResult$Results$Interval,binnedDataResult$Results$Treatment,binnedDataResult$Results$MeanDuration)
    names(analysis)<-c("Interval","Treatment","Y")
  }
  else if(Type=="MinInt") {
    gp<-ggplot(tmp,aes(x=Minutes,y=MinInt,color=Treatment,group=Treatment)) + 
      geom_errorbar(aes(ymin=MinInt-MinIntSEM, ymax=MinInt+MinIntSEM,color=Treatment), width=.1, position=pd) +
      geom_line(position=pd,size=1) +
      geom_point(position=pd, size=4, shape=21, fill="white") +xlab("Minutes") + ylab("Min Event Intensity")
    filename<-paste("BinnedMinIntPlots.pdf",sep="")
    analysis<-data.frame(binnedDataResult$Results$Interval,binnedDataResult$Results$Treatment,binnedDataResult$Results$MinInt)
    names(analysis)<-c("Interval","Treatment","Y")
  }
  else if(Type=="TimeBtw") {
    gp<-ggplot(tmp,aes(x=Minutes,y=MeanTimeBtw,color=Treatment,group=Treatment)) + 
      geom_errorbar(aes(ymin=MeanTimeBtw-MeanTimeBtwSEM, ymax=MeanTimeBtw+MeanTimeBtwSEM,color=Treatment), width=.1, position=pd) +
      geom_line(position=pd,size=1) +
      geom_point(position=pd, size=4, shape=21, fill="white") +xlab("Minutes") + ylab("Time Between Events (sec)")
    filename<-paste("BinnedMeanTimeBtwPlots.pdf",sep="")
    analysis<-data.frame(binnedDataResult$Results$Interval,binnedDataResult$Results$Treatment,binnedDataResult$Results$MeanTimeBtw)
    names(analysis)<-c("Interval","Treatment","Y")
  }
  else if(Type=="PI") {
    stop("Plot type not supported for single-well chambers.")
  }
  else if(Type=="EventPI") {
    stop("Plot type not supported for single-well chambers.")
  }
  else {
    stop("Plot type does not exist.")
  }
  show(gp)
  if(SaveToFile==TRUE){
    ggsave(filename,gp)
  }
  
  l<-lapply(split(analysis, analysis$Interval), aov, formula=Y ~ Treatment)
  cat("** Interval specific ANOVA results **\n\n")
  lapply(l,summary)
  
}
BinnedPlot.TwoWell.Trt<-function(binnedDataResult,Type="Licks",SaveToFile=FALSE){
  stats<-binnedDataResult$Stats
  results<-binnedDataResult$Results
  
  tmp<-c("LicksA","EventsA","MeanDurationA","MedDurationA","MeanTimeBtwA","MedTimeBtwA","MeanIntA","MedianIntA","MinIntA","MaxIntA")
  tmp2<-paste(tmp,"SEM",sep="")
  a.wells<-c("Treatment","Minutes",tmp,tmp2)
  
  tmp<-c("LicksB","EventsB","MeanDurationB","MedDurationB","MeanTimeBtwB","MedTimeBtwB","MeanIntB","MedianIntB","MinIntB","MaxIntB")
  tmp2<-paste(tmp,"SEM",sep="")
  b.wells<-c("Treatment","Minutes",tmp,tmp2)
  
  tmp<-c("Licks","Events","MeanDuration","MedDuration","MeanTimeBtw","MedTimeBtw","MeanInt","MedianInt","MinInt","MaxInt")
  tmp2<-paste(tmp,"SEM",sep="")
  new.names<-c("Treatment","Minutes",tmp,tmp2,"Well")
  tmpA<-data.frame(stats[,a.wells],rep("WellA",nrow(stats)))
  tmpB<-data.frame(stats[,b.wells],rep("WellB",nrow(stats)))
  names(tmpA)<-new.names
  names(tmpB)<-new.names
  newData<-rbind(tmpA,tmpB)
  
  pd <- position_dodge(5) # move them .05 to the left and right
  
  if(Type=="Licks") {
    ylabel<-"Licks"
    gp<-ggplot(newData,aes(x=Minutes,y=Licks,color=Treatment,group=Treatment)) + 
      geom_errorbar(aes(ymin=Licks-LicksSEM, ymax=Licks+LicksSEM,color=Treatment), width=.1, position=pd) +
      geom_line(position=pd,size=1) + facet_wrap(~Well)+
      geom_point(position=pd, size=4, shape=21, fill="white") +xlab("Minutes") + ylab(ylabel)
    filename<-paste("BinnedLicksPlots.pdf",sep="")
    analysis<-data.frame(results$Interval,results$Treatment,results$LicksA,results$LicksB)
    names(analysis)<-c("Interval","Treatment","YA","YB")
  }
  else if(Type=="Events") {
    ylabel<-"Events"
    gp<-ggplot(newData,aes(x=Minutes,y=Events,color=Treatment,group=Treatment)) + 
      geom_errorbar(aes(ymin=Events-EventsSEM, ymax=Events+EventsSEM,color=Treatment), width=.1, position=pd) +
      geom_line(position=pd,size=1) + facet_wrap(~Well)+
      geom_point(position=pd, size=4, shape=21, fill="white") +xlab("Minutes") + ylab("Events")
    filename<-paste("BinnedEventsPlots.pdf",sep="")
    analysis<-data.frame(results$Interval,results$Treatment,results$EventsA,results$EventsB)
    names(analysis)<-c("Interval","Treatment","YA","YB")
  }
  else if(Type=="Durations") {
    gp<-ggplot(newData,aes(x=Minutes,y=MeanDuration,color=Treatment,group=Treatment)) + 
      geom_errorbar(aes(ymin=MeanDuration-MeanDurationSEM, ymax=MeanDuration+MeanDurationSEM,color=Treatment), width=.1, position=pd) +
      geom_line(position=pd,size=1) + facet_wrap(~Well)+
      geom_point(position=pd, size=4, shape=21, fill="white") +xlab("Minutes") + ylab("Event Duration (sec)")
    filename<-paste("BinnedDurationPlots.pdf",sep="")
    analysis<-data.frame(results$Interval,results$Treatment,results$MeanDurationA,results$MeanDurationB)
    names(analysis)<-c("Interval","Treatment","YA","YB")
  }
  else if(Type=="MinInt") {
    gp<-ggplot(newData,aes(x=Minutes,y=MinInt,color=Treatment,group=Treatment)) + 
      geom_errorbar(aes(ymin=MinInt-MinIntSEM, ymax=MinInt+MinIntSEM,color=Treatment), width=.1, position=pd) +
      geom_line(position=pd,size=1) + facet_wrap(~Well)+
      geom_point(position=pd, size=4, shape=21, fill="white") +xlab("Minutes") + ylab("Min Event Intensity")
    filename<-paste("BinnedMinIntPlots.pdf",sep="")
    analysis<-data.frame(results$Interval,results$Treatment,results$MinIntA,results$MinIntB)
    names(analysis)<-c("Interval","Treatment","YA","YB")
  }
  else if(Type=="TimeBtw") {
    gp<-ggplot(newData,aes(x=Minutes,y=MeanTimeBtw,color=Treatment,group=Treatment)) + 
      geom_errorbar(aes(ymin=MeanTimeBtw-MeanTimeBtwSEM, ymax=MeanTimeBtw+MeanTimeBtwSEM,color=Treatment), width=.1, position=pd) +
      geom_line(position=pd,size=1) + facet_wrap(~Well)+
      geom_point(position=pd, size=4, shape=21, fill="white") +xlab("Minutes") + ylab("Time Between Events (sec)")
    filename<-paste("BinnedMeanTimeBtwPlots.pdf",sep="")
    analysis<-data.frame(results$Interval,results$Treatment,results$MeanTimeBtwA,results$MeanTimeBtwB)
    names(analysis)<-c("Interval","Treatment","YA","YB")
  }
  else if(Type=="PI"){
    Licks<-stats$LicksA+stats$LicksB
    stats<-data.frame(stats,Licks)
    gp<-ggplot(stats,aes(x=Minutes,y=PI,color=Treatment,group=Treatment)) + 
      geom_errorbar(aes(ymin=PI-PISEM, ymax=PI+PISEM,color=Treatment), width=.1, position=pd) +
      geom_line(position=pd,size=1) + ylim(c(-1.05,1.05)) +
      geom_point(position=pd, size=4, shape=21, aes(fill = Licks)) +xlab("Minutes") + ylab("PI (Licks)")
    filename<-paste("BinnedPIPlots.pdf",sep="")
    analysis<-data.frame(results$Interval,results$Treatment,results$PI,results$PI)
    names(analysis)<-c("Interval","Treatment","YA","YB")
  }
  else if(Type=="EventPI"){
    Events<-stats$EventsA+stats$EventsB
    stats<-data.frame(stats,Events)
    gp<-ggplot(stats,aes(x=Minutes,y=EventPI,color=Treatment,group=Treatment)) + 
      geom_errorbar(aes(ymin=EventPI-EventPISEM, ymax=EventPI+EventPISEM,color=Treatment), width=.1, position=pd) +
      geom_line(position=pd,size=1) + ylim(c(-1.05,1.05)) +
      geom_point(position=pd, size=4, shape=21,  aes(fill = Events)) +xlab("Minutes") + ylab("Event PI")
    filename<-paste("BinnedEventPIPlots.pdf",sep="")
    analysis<-data.frame(results$Interval,results$Treatment,results$EventPI,results$EventPI)
    names(analysis)<-c("Interval","Treatment","YA","YB")
  }
  else {
    stop("Plot type does not exist.")
  }
  show(gp)
  if(SaveToFile==TRUE){
    ggsave(filename,gp)
  }
  l<-lapply(split(analysis, analysis$Interval), aov, formula=YA ~ Treatment)
  cat("\n\n\n** Interval specific ANOVA results for Well A **\n\n")
  print(lapply(l,summary))
  
  l<-lapply(split(analysis, analysis$Interval), aov, formula=YB ~ Treatment)
  cat("\n\n\n** Interval specific ANOVA results for Well B **\n\n")
  print(lapply(l,summary))
}
SimpleDataPlot.OneWell<-function(summaryResults,Type="Licks",SaveToFile=FALSE){
  results<-summaryResults$Results
  results<-subset(results,Treatment!="None")
  
  if(Type=="Licks") {
    filename<-paste("SimpleLicksPlot.pdf",sep="")
    ylabel<-"Licks"
    r<-"Licks"
    gp<-(ggplot(results, aes(Treatment, Licks)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
           ylim(c(min(results$Licks),max(results$Licks))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE))
    analysis<-data.frame(results$Treatment,results$Licks)
    names(analysis)<-c("Treatment","Y")
  }
  else if(Type=="Events") {
    filename<-paste("SimpleEventsPlot.pdf",sep="")
    ylabel<-"Events"
    r<-"Events"
    gp<-(ggplot(results, aes(Treatment, Events)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
           ylim(c(min(results$Events),max(results$Events))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE))
    analysis<-data.frame(results$Treatment,results$Events)
    names(analysis)<-c("Treatment","Y")
  }
  else if(Type=="Durations") {
    filename<-paste("SimpleDurationsPlot.pdf",sep="")
    ylabel<-"Duration (sec)"
    r<-"Duration"
    gp<-(ggplot(results, aes(Treatment, MeanDuration)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
           ylim(c(min(results$MeanDuration),max(results$MeanDuration))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE))
    analysis<-data.frame(results$Treatment,results$MeanDuration)
    names(analysis)<-c("Treatment","Y")
    
  }
  else if(Type=="MinInt") {
    filename<-paste("SimpleMinIntPlot.pdf",sep="")
    ylabel<-"Minimum Event Intensity"
    r<-"Min Intensity"
    gp<-(ggplot(results, aes(Treatment, MinInt)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
           ylim(c(min(results$MinInt),max(results$MinInt))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE))
    analysis<-data.frame(results$Treatment,results$MinInt)
    names(analysis)<-c("Treatment","Y")
    
  }
  else if(Type=="TimeBtw") {
    filename<-paste("SimpleMeantTimeBtwPlot.pdf",sep="")
    ylabel<-"Time Between Events (sec)"
    r<-"Time Btw"
    gp<-(ggplot(results, aes(Treatment, MeanTimeBtw)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
           ylim(c(min(results$MeanTimeBtw),max(results$MeanTimeBtw))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE))
    analysis<-data.frame(results$Treatment,results$MeanTimeBtw)
    names(analysis)<-c("Treatment","Y")
  }
  else if(Type=="PI") {
    stop("Plot type not supported for single-well chambers.")
  }
  else if(Type=="EventPI") {
    stop("Plot type not supported for single-well chambers.")
  }
  else {
    stop("Plot type does not exist.")
  }
  show(gp)
  if(SaveToFile==TRUE)
    ggsave(filename = filename)
  r2<-paste("\n** ",r," **\n")
  cat(r2)
  print(summary(aov(Y~Treatment,data=analysis)))
}
SimpleDataPlot.TwoWell<-function(summaryResults,Type="Licks",SaveToFile=FALSE){
  results<-summaryResults$Results
  results<-subset(results,Treatment!="None")
  
  a.wells<-c("Treatment","LicksA","EventsA","MeanDurationA","MedDurationA","MeanTimeBtwA","MedTimeBtwA","MeanIntA","MedianIntA","MinIntA","MaxIntA")
  b.wells<-c("Treatment","LicksB","EventsB","MeanDurationB","MedDurationB","MeanTimeBtwB","MedTimeBtwB","MeanIntB","MedianIntB","MinIntB","MaxIntB")
  new.names<-c("Treatment","Licks","Events","MeanDuration","MedDuration","MeanTimeBtw","MedTimeBtw","MeanInt","MedianInt","MinInt","MaxInt","Well")
  tmpA<-data.frame(results[,a.wells],rep("WellA",nrow(results)))
  tmpB<-data.frame(results[,b.wells],rep("WellB",nrow(results)))
  names(tmpA)<-new.names
  names(tmpB)<-new.names
  newData<-rbind(tmpA,tmpB)
  
  if(Type=="Licks") {
    filename<-paste("SimpleLicksPlot.pdf",sep="")
    ylabel<-"Licks"
    r<-"Licks"
    gp<-ggplot(newData, aes(Treatment, Licks)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) + facet_wrap(~Well)+
      ylim(c(min(newData$Licks),max(newData$Licks))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE)
    analysis<-data.frame(results$Treatment,results$LicksA,results$LicksB)
    names(analysis)<-c("Treatment","YA","YB")
  }
  else if(Type=="Events") {
    filename<-paste("SimpleEventsPlot.pdf",sep="")
    ylabel<-"Events"
    r<-"Events"
    gp<-ggplot(newData, aes(Treatment, Events)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) + facet_wrap(~Well)+
      ylim(c(min(newData$Events),max(newData$Events))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE)
    analysis<-data.frame(results$Treatment,results$EventsA,results$EventsB)
    names(analysis)<-c("Treatment","YA","YB")
  }
  else if(Type=="Durations") {
    filename<-paste("SimpleDurationsPlot.pdf",sep="")
    ylabel<-"Duration (sec))"
    r<-"Duration"
    gp<-ggplot(newData, aes(Treatment, MeanDuration)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) + facet_wrap(~Well)+
      ylim(c(min(newData$MeanDuration),max(newData$MeanDuration))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE)
    analysis<-data.frame(results$Treatment,results$MeanDurationA,results$MeanDurationB)
    names(analysis)<-c("Treatment","YA","YB")
    
  }
  else if(Type=="MinInt") {
    filename<-paste("SimpleMinIntPlot.pdf",sep="")
    ylabel<-"Minimum Event Intensity"
    r<-"Min Intensity"
    gp<-ggplot(newData, aes(Treatment,MinInt)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) + facet_wrap(~Well)+
      ylim(c(min(newData$MinInt),max(newData$MinInt))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE)
    analysis<-data.frame(results$Treatment,results$MinIntA,results$MinIntB)
    names(analysis)<-c("Treatment","YA","YB")
  }
  else if(Type=="TimeBtw") {
    filename<-paste("SimpleMeantTimeBtwPlot.pdf",sep="")
    ylabel<-"Time Between Events (sec) (wmean(A,B))"
    r<-"Time Btw"
    gp<-ggplot(newData, aes(Treatment, MeanTimeBtw)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) + facet_wrap(~Well)+
      ylim(c(min(newData$MeanTimeBtw),max(newData$MeanTimeBtw))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE)
    analysis<-data.frame(results$Treatment,results$MeanTimeBtwA,results$MeanTimeBtwB)
    names(analysis)<-c("Treatment","YA","YB")
  }
  else if(Type=="PI") {
    filename<-paste("SimplePIPlot.pdf",sep="")
    ylabel<-"Licks PI"
    
    Licks<-results$LicksA+results$LicksB
    results<-data.frame(results,Licks)
    results<-results[Licks>0,]
    
    r<-"PI (Licks)"
    gp<-ggplot(results, aes(results$Treatment, results$PI)) + geom_boxplot(aes(fill = results$Treatment),outlier.size=-1) + geom_jitter(aes(color = Licks),size=3,height=0) +
      ylim(c(-1.05,1.05)) + ggtitle(r) + xlab("Treatment") +ylab("PI") + guides(fill=FALSE)
    
    analysis<-data.frame(results$Treatment,results$PI,results$PI)
    names(analysis)<-c("Treatment","YA","YB")
  }
  else if(Type=="EventPI") {
    filename<-paste("SimpleEventPIPlot.pdf",sep="")
    ylabel<-"Event PI"
    
    Events<-results$EventsA+results$EventsB
    results<-data.frame(results,Events)
    results<-results[Events>0,]
    
    r<-"PI (Events)"
    gp<-ggplot(results, aes(results$Treatment, results$EventPI)) + geom_boxplot(aes(fill = results$Treatment),outlier.size=-1) + geom_jitter(aes(color = Events),size=3,height=0) +
      ylim(c(-1.05,1.05)) + ggtitle(r) + xlab("Treatment") +ylab("Event PI") + guides(fill=FALSE)
    
    analysis<-data.frame(results$Treatment,results$EventPI,results$EventPI)
    names(analysis)<-c("Treatment","YA","YB")
  }
  else {
    stop("Plot type does not exist.")
  }
  show(gp)
  if(SaveToFile==TRUE)
    ggsave(filename = filename)
  
  r2<-paste("\n** ",r," **\n")
  cat(r2)
  cat("\nANOVA results for Well A\n")
  print(summary(aov(YA~Treatment,data=analysis)))
  
  cat("\n\nANOVA results for Well B\n")
  print(summary(aov(YB~Treatment,data=analysis)))
}

## Functions used for cumulative division plots of basic stats.
DivisionPlots.OneWell<-function(monitors,parameters,expDesign,range=c(0,0),divisions=1,Type="Licks",SaveToFile=FALSE,TransformLicks=TRUE){
  ranges<-matrix(rep(NA,divisions*2),ncol=2)
  if(divisions==1)
    ranges[1,]<-range
  else {
    if(range[2]==0){
      ## Come back to here
      dfm<-GetDFM(monitors[1])
      if(sum(is.na(dfm)))
        stop("DFM Missing")
      last.time<-LastSampleData(dfm)$Minutes
      breaks<-seq(from=range[1], to=last.time, length=divisions+1)
      ranges.1<-range[1]
      ranges.2<-breaks[-1]
      ranges<-round(cbind(ranges.1,ranges.2))
    }
    else {
      breaks<-seq(from=range[1], to=range[2], length=divisions+1)
      ranges.1<-range[1]
      ranges.2<-breaks[-1]
      ranges<-round(cbind(ranges.1,ranges.2))
    }
  }
  if(SaveToFile==TRUE){
    filename<-paste("DivisionPlots_",Type,".pdf",sep="")
    pdf(file=filename)
  }
  if(divisions==1) {
    tmp<-Feeding.Summary.Monitors(monitors,parameters,expDesign,range,FALSE,TransformLicks)
    results<-tmp$Results
    results<-subset(results,Treatment!="None")
    
    if(Type=="Licks"){
      if(TransformLicks==TRUE){
        r<-paste("Transformed Licks -- Range(min): (",range[1],",",range[2],")",sep="")
        ylabel<-"Transformed Licks"
      }
      else {
        r<-paste("Licks -- Range(min): (",range[1],",",range[2],")",sep="")
        ylabel<-"Licks"
      }
      
      print(ggplot(results, aes(Treatment, Licks)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
              ylim(c(min(results$Licks),max(results$Licks))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE))
      r2<-paste("\n** ",r," **\n")
      cat(r2)
      print(summary(aov(Licks~Treatment,data=results)))
    }
    else if(Type=="Events"){
      r<-paste("Events-Range: (",range[1],",",range[2],")",sep="")
      print(ggplot(results, aes(Treatment, Events)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
              ylim(c(min(results$Events),max(results$Events))) + ggtitle(r) + xlab("Treatment") +ylab("Events") + guides(fill=FALSE))
      cat("** ANOVA results Full Range **\n\n")
      print(summary(aov(Events~Treatment,data=results)))
    }
    else if(Type=="Durations"){
      r<-paste("Durations -- Range(min): (",range[1],",",range[2],")",sep="")
      print(ggplot(results, aes(Treatment, MeanDuration)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
              ylim(c(min(results$MeanDuration),max(results$MeanDuration))) + ggtitle(r) + xlab("Treatment") +ylab("Mean Duration (sec)") + guides(fill=FALSE))
      print(summary(aov(MeanDuration~Treatment,data=results)))
    }
    else if(Type=="TimeBtw"){
      r<-paste("Time Btw -- Range(min): (",range[1],",",range[2],")",sep="")
      print(ggplot(results, aes(Treatment, MeanTimeBtw)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
              ylim(c(min(results$MeanTimeBtw),max(results$MeanTimeBtw))) + ggtitle(r) + xlab("Treatment") +ylab("Mean Time Between Events (sec)") + guides(fill=FALSE))
      print(summary(aov(MeanTimeBtw~Treatment,data=results)))
    }
    else if(Type=="MinInt"){
      r<-paste("Min Intensity -- Range(min): (",range[1],",",range[2],")",sep="")
      print(ggplot(results, aes(Treatment, MinInt)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
              ylim(c(min(results$MinInt),max(results$MinInt))) + ggtitle(r) + xlab("Treatment") +ylab("Min Event Intensity") + guides(fill=FALSE))
      print(summary(aov(MinInt~Treatment,data=results)))
    }
    else if(Type=="PI"){
      stop("Plot type not supported for single-well chambers.")
    }
    else if(Type=="EventPI"){
      stop("Plot type not supported for single-well chambers.")
    }
    else {
      if(SaveToFile==TRUE){
        graphics.off()
      }
      stop("Plot type does not exist.")
    }
  }
  else {
    p<-list()
    for(i in 1:divisions)
      local({
        tmp<-Feeding.Summary.Monitors(monitors,parameters,expDesign,ranges[i,],FALSE,TransformLicks)
        results<-tmp$Results
        results<-subset(results,Treatment!="None")
        if(Type=="Licks"){
          if(TransformLicks==TRUE){
            r<-paste("Transformed Licks -- Range(min): (",ranges[i,1],",",ranges[i,2],")",sep="")
            ylabel<-"Transformed Licks"
          }
          else {
            r<-paste("Licks -- Range(min): (",ranges[i,1],",",ranges[i,2],")",sep="")
            ylabel<-"Licks"
          }
          r2<-paste("\n** ",r," **\n")
          cat(r2)
          print(summary(aov(Licks~Treatment,data=results)))
          p[[i]]<<-(ggplot(results, aes(Treatment, Licks)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
                      ylim(c(min(results$Licks),max(results$Licks))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE))
        }
        else if(Type=="Events"){
          r<-paste("Events-Range: (",ranges[i,1],",",ranges[i,2],")",sep="")
          p[[i]]<<-(ggplot(results, aes(Treatment, Events)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
                      ylim(c(min(results$Events),max(results$Events))) + ggtitle(r) + xlab("Treatment") +ylab("Events") + guides(fill=FALSE))
          ppp<-paste("\n** ANOVA results -- Range(min): (",ranges[i,1],",",ranges[i,2],")",sep="")
          cat(ppp)
          print(summary(aov(Events~Treatment,data=results)))
        }
        else if(Type=="Durations"){
          r<-paste("Durations -- Range(min): (",ranges[i,1],",",ranges[i,2],")",sep="")
          p[[i]]<<-(ggplot(results, aes(Treatment, MeanDuration)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
                      ylim(c(min(results$MeanDuration),max(results$MeanDuration))) + ggtitle(r) + xlab("Treatment") +ylab("Mean Duration (sec)") + guides(fill=FALSE))
          print(summary(aov(MeanDuration~Treatment,data=results)))
        }
        else if(Type=="TimeBtw"){
          r<-paste("Time Btw -- Range(min): (",ranges[i,1],",",ranges[i,2],")",sep="")
          p[[i]]<<-(ggplot(results, aes(Treatment, MeanTimeBtw)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
                      ylim(c(min(results$MeanTimeBtw),max(results$MeanTimeBtw))) + ggtitle(r) + xlab("Treatment") +ylab("Mean Time Between Events (sec)") + guides(fill=FALSE))
          print(summary(aov(MeanTimeBtw~Treatment,data=results)))
        }
        else if(Type=="MinInt"){
          r<-paste("Min Intensity -- Range(min): (",ranges[i,1],",",ranges[i,2],")",sep="")
          p[[i]]<<-(ggplot(results, aes(Treatment, MinInt)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
                      ylim(c(min(results$MinInt),max(results$MinInt))) + ggtitle(r) + xlab("Treatment") +ylab("Min Event Intensity") + guides(fill=FALSE))
          print(summary(aov(MinInt~Treatment,data=results)))
        }
        else if(Type=="PI"){
          stop("Plot type not supported for single-well chambers.")
        }
        else if(Type=="EventPI"){
          stop("Plot type not supported for single-well chambers.")
        }
        else {
          if(SaveToFile==TRUE){
            graphics.off()
          }
          stop("Plot type does not exist.")
        }
      })
    if(divisions<5)
      numcols<-2
    else if(divisions<10)
      numcols<-3
    else if(divisions<17)
      numcols<-4
    else
      numcols<-5
    multiplot(plotlist=p,cols=numcols)
  }
  if(SaveToFile==TRUE)
    graphics.off()
}
DivisionPlots.TwoWell<-function(monitors,parameters,expDesign,range=c(0,0),divisions=1,Type="Licks",SaveToFile=FALSE,TransformLicks=TRUE){
  ranges<-matrix(rep(NA,divisions*2),ncol=2)
  if(divisions==1)
    ranges[1,]<-range
  else {
    if(range[2]==0){
      dfm<-GetDFM(monitors[1])
      if(sum(is.na(dfm)))
        stop("DFM Missing")
      last.time<-LastSampleData(dfm)$Minutes
      breaks<-seq(from=range[1], to=last.time, length=divisions+1)
      ranges.1<-range[1]
      ranges.2<-breaks[-1]
      ranges<-round(cbind(ranges.1,ranges.2))
    }
    else {
      breaks<-seq(from=range[1], to=range[2], length=divisions+1)
      ranges.1<-range[1]
      ranges.2<-breaks[-1]
      ranges<-round(cbind(ranges.1,ranges.2))
    }
  }
  
  if(SaveToFile==TRUE){
    filename<-paste("DivisionPlots_",Type,".pdf",sep="")
    pdf(file=filename)
  }
  
  if(divisions==1) {
    tmp<-Feeding.Summary.Monitors(monitors,parameters,expDesign,range,FALSE,TransformLicks)
    results<-tmp$Results
    results<-subset(results,Treatment!="None")
    
    if(Type=="Licks"){
      Licks<-results$LicksA+results$LicksB
      if(TransformLicks==TRUE){
        r<-paste("Transformed Total Licks -- Range(min): (",range[1],",",range[2],")",sep="")
        ylabel<-"Transformed Total Licks (A+B)"
      }
      else {
        r<-paste("Licks -- Range(min): (",range[1],",",range[2],")",sep="")
        ylabel<-"Total Licks (A+B)"
      }
      results<-data.frame(results,Licks)
      r2<-paste("\n** ",r," **\n")
      cat(r2)
      print(summary(aov(Licks~Treatment,data=results)))
      print(ggplot(results, aes(Treatment, Licks)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
              ylim(c(min(results$Licks),max(results$Licks))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE))
    }
    else if(Type=="Events"){
      Events<-results$EventsA+results$EventsB
      results<-data.frame(results,Events)
      r<-paste("Events-Range: (",range[1],",",range[2],")",sep="")
      print(ggplot(results, aes(Treatment, Events)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
              ylim(c(min(results$Events),max(results$Events))) + ggtitle(r) + xlab("Treatment") +ylab("Total Events (A+B)") + guides(fill=FALSE))
      cat("** ANOVA results Full Range **\n\n")
      print(summary(aov(Events~Treatment,data=results)))
    }
    else if(Type=="Durations"){
      MeanDuration <- ((results$MeanDurationA*results$EventsA)+ (results$MeanDurationB*results$EventsB))/(results$EventsA+results$EventsB)
      results<-data.frame(results,MeanDuration)
      r<-paste("Durations -- Range(min): (",range[1],",",range[2],")",sep="")
      print(ggplot(results, aes(Treatment, MeanDuration)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
              ylim(c(min(results$MeanDuration),max(results$MeanDuration))) + ggtitle(r) + xlab("Treatment") +ylab("Mean Duration (A and B)") + guides(fill=FALSE))
      print(summary(aov(MeanDuration~Treatment,data=results)))
    }
    else if(Type=="TimeBtw"){
      MeanTimeBtw <- ((results$MeanTimeBtwA*results$EventsA)+ (results$MeanTimeBtwB*results$EventsB))/(results$EventsA+results$EventsB)
      results<-data.frame(results,MeanTimeBtw)
      r<-paste("Time Btw -- Range(min): (",range[1],",",range[2],")",sep="")
      print(ggplot(results, aes(Treatment, MeanTimeBtw)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
              ylim(c(min(results$MeanTimeBtw),max(results$MeanTimeBtw))) + ggtitle(r) + xlab("Treatment") +ylab("Time Between") + guides(fill=FALSE))
      print(summary(aov(MeanTimeBtw~Treatment,data=results)))
    }
    else if(Type=="MinInt"){
      MinInt<-pmin(results$MinIntA,results$MinIntB)
      results<-data.frame(results,MinInt)
      r<-paste("Min Intensity -- Range(min): (",range[1],",",range[2],")",sep="")
      print(ggplot(results, aes(Treatment, MinInt)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
              ylim(c(min(results$MinInt),max(results$MinInt))) + ggtitle(r) + xlab("Treatment") +ylab("Min Event Intensity") + guides(fill=FALSE))
      print(summary(aov(MinInt~Treatment,data=results)))
    }
    else if(Type=="PI"){
      Licks<-results$LicksA+results$LicksB
      results<-data.frame(results,Licks)
      results<-results[Licks>0,]
      r<-paste("PI -- Range: (",range[1],",",range[2],")",sep="")
      print(ggplot(results, aes(Treatment, PI)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(aes(color = Licks),size=3,height=0) +
              ylim(c(-1.05,1.05)) + ggtitle(r) + xlab("Treatment") +ylab("PI") + guides(fill=FALSE))
      print(summary(aov(PI~Treatment,data=results)))
    }
    else if(Type=="EventPI"){
      Events<-results$EventsA+results$EventsB
      results<-data.frame(results,Events)
      results<-results[Events>0,]
      r<-paste("Event PI -- Range: (",range[1],",",range[2],")",sep="")
      print(ggplot(results, aes(Treatment, EventPI)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(aes(color = Events),size=3,height=0) +
              ylim(c(-1.05,1.05)) + ggtitle(r) + xlab("Treatment") +ylab("Event PI") + guides(fill=FALSE))
      print(summary(aov(EventPI~Treatment,data=results)))
    }
    else {
      if(SaveToFile==TRUE){
        graphics.off()
      }
      stop("Plot type does not exist.")
    }
  }
  else {
    p<-list()
    for(i in 1:divisions)
      local({
        ## Because we are adding licks together, we should transform only after adding!!
        ## So this function should always work on non-transformed data.
        tmp<-Feeding.Summary.Monitors(monitors,parameters,expDesign,ranges[i,],FALSE,TransformLicks=FALSE)
        results<-tmp$Results
        results<-subset(results,Treatment!="None")
        if(Type=="Licks"){
          Licks<-results$LicksA+results$LicksB
          if(TransformLicks==TRUE){
            r<-paste("Transformed Total Licks -- Range(min): (",ranges[i,1],",",ranges[i,2],")"
                     ,sep="")
            ylabel<-"Transformed Totoal Licks (A+B)"
            Licks<-Licks^0.25
          }
          else {
            r<-paste("Total Licks -- Range(min): (",ranges[i,1],",",ranges[i,2],")",sep="")
            ylabel<-"Total Licks (A+B)"
          }
          results<-data.frame(results,Licks)
          r2<-paste("\n** ",r," **\n")
          cat(r2)
          print(summary(aov(Licks~Treatment,data=results)))
          p[[i]]<<-(ggplot(results, aes(Treatment, Licks)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
                      ylim(c(min(results$Licks),max(results$Licks))) + ggtitle(r) + xlab("Treatment") +ylab(ylabel) + guides(fill=FALSE))
        }
        else if(Type=="Events"){
          Events<-results$EventsA+results$EventsB
          results<-data.frame(results,Events)
          r<-paste("Events-Range: (",ranges[i,1],",",ranges[i,2],")",sep="")
          p[[i]]<<-(ggplot(results, aes(Treatment, Events)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
                      ylim(c(min(results$Events),max(results$Events))) + ggtitle(r) + xlab("Treatment") +ylab("Total Events (A+B)") + guides(fill=FALSE))
          ppp<-paste("\n** ANOVA results -- Range(min): (",ranges[i,1],",",ranges[i,2],")",sep="")
          cat(ppp)
          print(summary(aov(Events~Treatment,data=results)))
        }
        else if(Type=="Durations"){
          MeanDuration <- ((results$MeanDurationA*results$EventsA)+ (results$MeanDurationB*results$EventsB))/(results$EventsA+results$EventsB)
          results<-data.frame(results,MeanDuration)
          r<-paste("Durations -- Range(min): (",ranges[i,1],",",ranges[i,2],")",sep="")
          p[[i]]<<-(ggplot(results, aes(Treatment, MeanDuration)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
                      ylim(c(min(results$MeanDuration),max(results$MeanDuration))) + ggtitle(r) + xlab("Treatment") +ylab("Mean Duration (A and B") + guides(fill=FALSE))
          print(summary(aov(MeanDuration~Treatment,data=results)))
        }
        else if(Type=="TimeBtw"){
          MeanTimeBtw <- ((results$MeanTimeBtwA*results$EventsA)+ (results$MeanTimeBtwB*results$EventsB))/(results$EventsA+results$EventsB)
          results<-data.frame(results,MeanTimeBtw)
          r<-paste("Time Btw -- Range(min): (",ranges[i,1],",",ranges[i,2],")",sep="")
          p[[i]]<<-(ggplot(results, aes(Treatment,MeanTimeBtw)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
                      ylim(c(min(results$MeanTimeBtw),max(results$MeanTimeBtw))) + ggtitle(r) + xlab("Treatment") +ylab("Time Between") + guides(fill=FALSE))
          print(summary(aov(MeanTimeBtw~Treatment,data=results)))
        }
        else if(Type=="MinInt"){
          MinInt<-pmin(results$MinIntA,results$MinIntB)
          results<-data.frame(results,MinInt)
          r<-paste("Min Intensity -- Range(min): (",ranges[i,1],",",ranges[i,2],")",sep="")
          p[[i]]<<-(ggplot(results, aes(Treatment, MinInt)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(size=3,height=0) +
                      ylim(c(min(results$MinInt),max(results$MinInt))) + ggtitle(r) + xlab("Treatment") +ylab("Min Event Intensity") + guides(fill=FALSE))
          print(summary(aov(MinInt~Treatment,data=results)))
        }
        else if(Type=="PI"){
          Licks<-results$LicksA+results$LicksB
          results<-data.frame(results,Licks)
          results<-results[Licks>0,]
          r<-paste("PI -- Range: (",ranges[i,1],",",ranges[i,2],")",sep="")
          p[[i]]<<-(ggplot(results, aes(Treatment, PI)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(aes(color = Licks),size=3,height=0) +
                      ylim(c(-1.05,1.05)) + ggtitle(r) + xlab("Treatment") +ylab("PI") + guides(fill=FALSE))
          print(summary(aov(PI~Treatment,data=results)))
        }
        else if(Type=="EventPI"){
          Events<-results$EventsA+results$EventsB
          results<-data.frame(results,Events)
          results<-results[Events>0,]
          r<-paste("Event PI -- Range: (",ranges[i,1],",",ranges[i,2],")",sep="")
          p[[i]]<<-(ggplot(results, aes(Treatment, EventPI)) + geom_boxplot(aes(fill = Treatment),outlier.size=-1) + geom_jitter(aes(color = Events),size=3,height=0) +
                      ylim(c(-1.05,1.05)) + ggtitle(r) + xlab("Treatment") +ylab("Event PI") + guides(fill=FALSE))
          print(summary(aov(EventPI~Treatment,data=results)))
        }
        else {
          if(SaveToFile==TRUE){
            graphics.off()
          }
          stop("Plot type does not exist.")
        }
      })
    if(divisions<5)
      numcols<-2
    else if(divisions<10)
      numcols<-3
    else if(divisions<17)
      numcols<-4
    else
      numcols<-5
    multiplot(plotlist=p,cols=numcols)
  }
  if(SaveToFile==TRUE)
    graphics.off()
}

## For the individual DFM plots.
PlotBins.Licks.DFM.OneWell<-function(dfm,binsize.min=30,range=c(0,0),TransformLicks=TRUE){
  if(dfm$Parameters$Chamber.Size!=1)
    stop("This function is for single chambers only")
  binnedData<-BinnedFeeding.Summary.DFM(dfm,binsize.min,range,TransformLicks)
  ylabel<-"Licks"
  xlabel<-"Minutes"
  ttl<-paste("DFM:",dfm$ID)
  if(TransformLicks==TRUE) {
    ylabel<-"Transformed Licks"
    ttl<-paste("DFM:",dfm$ID," (Transformed)")
  }
  gp<-ggplot(binnedData,aes(x=Minutes,y=Licks,fill=Chamber)) + geom_bar(stat="identity") + facet_grid(Chamber ~ .) + ggtitle(paste("DFM:",dfm$ID)) +
    theme(legend.position = "none") + ylab(ylabel) + xlab(xlabel)
  show(gp)
}
PlotBins.Licks.DFM.TwoWell<-function(dfm,binsize.min,range=c(0,0),TransformLicks=TRUE){
  if(dfm$Parameters$Chamber.Size!=2)
    stop("This function is for two-well chambers only")
  binnedData<-BinnedFeeding.Summary.DFM(dfm,binsize.min,range,TransformLicks)
  ylabel<-"Licks"
  xlabel<-"Minutes"
  ttl<-paste("DFM:",dfm$ID)
  if(TransformLicks==TRUE) {
    ylabel<-"Transformed Licks"
    ttl<-paste("DFM:",dfm$ID," (Transformed)")
  }
  tmp2<-melt(binnedData,id.vars=c("Minutes","Chamber"),measure.vars=c("LicksA","LicksB"))
  names(tmp2)[3]<-"Well"
  gp<-ggplot(tmp2,aes(x=Minutes,y=value,fill=Well)) + geom_bar(stat="identity") + facet_grid(Chamber ~ .) + ggtitle(paste("DFM:",dfm$ID)) +
    ylab(ylabel) + xlab(xlabel)
  show(gp)
}
PlotBins.Durations.DFM.OneWell<-function(dfm,binsize.min=30,range=c(0,0)){
  if(dfm$Parameters$Chamber.Size!=1)
    stop("This function is for single chambers only")
  binnedData<-BinnedFeeding.Summary.DFM(dfm,binsize.min,range,FALSE)
  ylabel<-"Avg Duration"
  xlabel<-"Minutes"
  ttl<-paste("DFM:",dfm$ID)
  gp<-ggplot(binnedData,aes(x=Minutes,y=Duration,fill=Chamber)) + geom_bar(stat="identity") + facet_grid(Chamber ~ .) + ggtitle(paste("DFM:",dfm$ID)) +
    theme(legend.position = "none") + ylab(ylabel) + xlab(xlabel)
  show(gp)
}
PlotBins.Durations.DFM.TwoWell<-function(dfm,binsize.min,range=c(0,0)){
  if(dfm$Parameters$Chamber.Size!=2)
    stop("This function is for two-well chambers only")
  binnedData<-BinnedFeeding.Summary.DFM(dfm,binsize.min,range,FALSE)
  ylabel<-"Avg Duration"
  xlabel<-"Minutes"
  ttl<-paste("DFM:",dfm$ID)
  tmp2<-melt(binnedData,id.vars=c("Minutes","Chamber"),measure.vars=c("DurationA","DurationB"))
  names(tmp2)[3]<-"Well"
  gp<-ggplot(tmp2,aes(x=Minutes,y=value,fill=Well)) + geom_bar(stat="identity") + facet_grid(Chamber ~ .) + ggtitle(paste("DFM:",dfm$ID)) +
    ylab(ylabel) + xlab(xlabel)
  show(gp)
}


##### Utility functions #####

## This function takes 2 vectors, one with the events
## above a minimal threshold (minvec) and one that
## specifies events that pass a more stringent threshold (maxvec).
## Contiguous events are only kept if at least one
## value in the event, which is defined by minvec, is above
## the higher threshold, which is defined by max vec
## z <- c(TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE)
## zz <- c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE)
## Get.Surviving.Events(z,zz) -> (2 0 0 0 0 0 3 0 0)
Get.Surviving.Events<-function(minvec,maxvec){
  tmp<-Get.Events(minvec)
  result<-tmp
  indices<-(1:length(minvec))[tmp>0]
  for(i in indices){
    tmp2<-maxvec[i:(i+(tmp[i]-1))]
    if(sum(tmp2)==0)
      result[i]<-0  
  }
  result
}

## This function is the reverse of Get.Events
## (2 0 0 0 1 0 3 0 0) -> c(TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE)
Expand.Events<-function(eventvec){
  result<-rep(FALSE,length(eventvec))
  indices<-(1:length(eventvec))[eventvec>0]
  for(i in indices){
    result[i:(i+eventvec[i]-1)]<-TRUE   
  }
  result
}

## These functions are helper functions for the basic calculations
# This function replaces continuing events with zero and make the first event of that
# episode equal to its duration.
## c(TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE) -> (2 0 0 0 1 0 3 0 0)
Get.Events<-function(z){
  tmp<-rle(z)
  result<-c(-1)
  for(i in 1:length(tmp$lengths)){
    if(tmp$values[i]){
      tmp2<-c(tmp$lengths[i],rep(0,tmp$lengths[i]-1))
      result<-c(result,tmp2)
    }
    else {
      tmp2<-c(rep(0,tmp$lengths[i]))
      result<-c(result,tmp2)
    }
  }
  result[-1]
}
Get.Events.And.Intensities<-function(z,data){
  z<-Get.Events(z)
  max.inten<-rep(0,length(z))
  min.inten<-rep(0,length(z))
  sum.inten<-rep(0,length(z))
  avg.inten<-rep(0,length(z))
  
  indices<-(1:length(z))[z>0]
  for(i in indices){
    tmp2<-data[i:(i+(z[i]-1))]
    max.inten[i]<-max(tmp2)
    min.inten[i]<-min(tmp2)
    sum.inten[i]<-sum(tmp2)
    avg.inten[i]<-mean(tmp2)    
  }
  result<-data.frame(z,min.inten,max.inten,sum.inten,avg.inten)
  names(result)<-c("FeedingEvent","MinIntensity","MaxIntensity","SumIntensity","MeanIntensity")
  result
}

## This function will take a TRUE/FALSE vector, assumed to be minthresholded lick data
## and it will "bridge" runs of FALSE of less than 'thresh' entries with TRUE
## e.g. c(TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE) -> c(TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,TRUE,TRUE)
## with thresh <= 3
Link.Events<-function(z,thresh){
  tmp<-rle(z)
  result<-c(FALSE)
  for(i in 1:length(tmp$lengths)){
    if(tmp$values[i]){
      tmp2<-rep(TRUE,tmp$lengths[i])
      result<-c(result,tmp2)
    }
    else {
      if(tmp$lengths[i]>thresh){
        tmp2<-rep(FALSE,tmp$lengths[i])
        result<-c(result,tmp2)
      }
      else {
        tmp2<-rep(TRUE,tmp$lengths[i])
        result<-c(result,tmp2)
      }
    }
  }
  result[-1]
}

# This function replaces continuing events with zero and make the first event of that
# episode equal to its duration.
## c(TRUE,TRUE,FALSE,FALSE,TRUE,FALSE,TRUE,TRUE,TRUE) -> (0 0 2 0 0 1 0 0 0)
Get.Intervals<-function(z){
  tmp<-rle(z)
  result<-c(-1)
  for(i in 1:length(tmp$lengths)){
    if(!tmp$values[i]){
      tmp2<-c(tmp$lengths[i],rep(0,tmp$lengths[i]-1))
      result<-c(result,tmp2)
    }
    else {
      tmp2<-c(rep(0,tmp$lengths[i]))
      result<-c(result,tmp2)
    }
  }
  result[-1]
}
mySEM<-function(x){
  if(is.factor(x)){
    tmp<-unique(as.character(x))
  }
  else{
    x<-x[!is.na(x)]
    tmp<-sqrt(var(x)/length(x))
  }
  tmp
}
GetTreatmentForChamber<-function(dfmNum,chamberNum,expdesign){
  tmp<-subset(expdesign,DFM==dfmNum)
  tmp<-subset(tmp,Chamber==chamberNum)
  if(nrow(tmp)!=1)
    return("None")
  else
    return(as.character(tmp$Treatment))
}

GetExpDesignForChamber<-function(dfmNum,chamberNum,expdesign){
  tmp<-subset(expdesign,DFM==dfmNum)
  tmp<-subset(tmp,Chamber==chamberNum)
  if(nrow(tmp)!=1)
    return(NA)
  else
    return(tmp)
}

gt<-function(df,expDesign){
  return("hi")
}
GetTreatmentForRow<-function(df,expdesign){
  tmp<-subset(expdesign,DFM==df["DFM"] & Chamber==df["Chamber"])
  if(nrow(tmp)!=1)
    return("None")
  else
    return(as.character(tmp$Treatment))
}
## Need to fix treatment assignments in single and choice experiments!!
AppendTreatmentonResultsFrame<-function(results,expdesign){
  isChamberThere<-"Chamber" %in% names(results)
  if(isChamberThere==FALSE)
    stop("Need a chamber to assign treatment.")
  Treatment<-rep(NA,nrow(results))
  for(i in 1:nrow(results)){
    Treatment[i]<-GetTreatmentForChamber(results$DFM[i],results$Chamber[i],expdesign)
  }
  n<-names(results)
  n<-c("Treatment",n)
  results<-cbind(Treatment,results)
  names(results)<-n
  results$Treatment<-factor(results$Treatment)
  results
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



###**********************************************************************************************************************
###* ParametersClass **************************************************************************************************
###* 
###* 

## First, the short parameters class
## It is instantiated with a set of initial values
ParametersClass=function(Baseline.Window.Minutes=3, Feeding.Threshold=20,Feeding.Minimum=10){
  Tasting.Interval=c(Feeding.Minimum,Feeding.Threshold)
  Feeding.Minevents=1
  Tasting.Minevents=1
  Samples.Per.Second=5
  Feeding.Event.Link.Gap=5
  Chamber.Sets=matrix(1:12,ncol=2,byrow=TRUE)
  Chamber.Size=2
  PI.Multiplier=1
  list(Baseline.Window.Minutes=Baseline.Window.Minutes,Feeding.Threshold=Feeding.Threshold,
       Feeding.Minimum=Feeding.Minimum,Tasting.Minimum=Tasting.Interval[1],
       Tasting.Maximum=Tasting.Interval[2],
       Feeding.Minevents=Feeding.Minevents,Tasting.Minevents=Tasting.Minevents,Samples.Per.Second=Samples.Per.Second,Chamber.Size=Chamber.Size,
       Chamber.Sets=Chamber.Sets,Feeding.Event.Link.Gap=Feeding.Event.Link.Gap,PI.Multiplier=PI.Multiplier)
}
ParametersClass.SingleWell=function(Baseline.Window.Minutes=3, Feeding.Threshold=20,Feeding.Minimum=10){
  Tasting.Interval=c(5, Feeding.Threshold)
  Feeding.Minevents=1
  Tasting.Minevents=1
  Samples.Per.Second=5
  Chamber.Sets=matrix(1:12,ncol=1,byrow=TRUE)
  Chamber.Size=1
  PI.Multiplier=0
  Feeding.Event.Link.Gap=5
  list(Baseline.Window.Minutes=Baseline.Window.Minutes,Feeding.Threshold=Feeding.Threshold,
       Feeding.Minimum=Feeding.Minimum,Tasting.Minimum=Tasting.Interval[1],
       Tasting.Maximum=Tasting.Interval[2],
       Feeding.Minevents=Feeding.Minevents,Tasting.Minevents=Tasting.Minevents,Samples.Per.Second=Samples.Per.Second,Chamber.Size=Chamber.Size,
       Chamber.Sets=Chamber.Sets,Feeding.Event.Link.Gap=Feeding.Event.Link.Gap,PI.Multiplier=PI.Multiplier)
}
ParametersClass.TwoWell=function(Baseline.Window.Minutes=3, Feeding.Threshold=20,Feeding.Minimum=10){
  Tasting.Interval=c(5,Feeding.Threshold)
  Feeding.Minevents=1
  Tasting.Minevents=1
  Samples.Per.Second=5
  Chamber.Sets=matrix(1:12,ncol=2,byrow=TRUE)
  Chamber.Size=2
  PI.Multiplier=1
  Feeding.Event.Link.Gap=5
  list(Baseline.Window.Minutes=Baseline.Window.Minutes,Feeding.Threshold=Feeding.Threshold,
       Feeding.Minimum=Feeding.Minimum,Tasting.Minimum=Tasting.Interval[1],
       Tasting.Maximum=Tasting.Interval[2],
       Feeding.Minevents=Feeding.Minevents,Tasting.Minevents=Tasting.Minevents,Samples.Per.Second=Samples.Per.Second,Chamber.Size=Chamber.Size,
       Chamber.Sets=Chamber.Sets,Feeding.Event.Link.Gap=Feeding.Event.Link.Gap,PI.Multiplier=PI.Multiplier)
}

## change the initial values using this function
SetParameter<-function(p,Baseline.Window.Minutes=NA,Feeding.Threshold=NA, Feeding.Minimum=NA, Tasting.Interval=c(NA,NA),
                       Feeding.Minevents=NA,Tasting.Minevents=NA,
                       Samples.Per.Sec=NA, Chamber.Size=NA, Feeding.Event.Link.Gap=NA,PI.Multiplier=NA){
  tmp.O<-options()
  options(warn=-1)
  ## Change only those that are listed
  if(!is.na(Baseline.Window.Minutes)) {
    p$Baseline.Window.Minutes=Baseline.Window.Minutes  
  }
  if(!is.na(Feeding.Threshold)) {
    p$Feeding.Threshold=Feeding.Threshold
  }
  if(!is.na(Feeding.Minimum)) {
    p$Feeding.Minimum=Feeding.Minimum
  }
  if(sum(is.na(Tasting.Interval))!=2) {
    p$Tasting.Minimum=Tasting.Interval[1]
    p$Tasting.Maximum=Tasting.Interval[2]
  }
  if(!is.na(Feeding.Minevents)){
    p$Feeding.Minevents=Feeding.Minevents
  }
  if(!is.na(Tasting.Minevents)){
    p$Tasting.Minevents=Tasting.Minevents
  }
  if(!is.na(Samples.Per.Sec)){
    p$Samples.Per.Second=Samples.Per.Sec
  }
  if(!is.na(Chamber.Size)){
    p$Chamber.Size=Chamber.Size
  }
  if(!is.na(PI.Multiplier)){
    p$PI.Multiplier=PI.Multiplier
  }
  if(!is.na(Feeding.Event.Link.Gap)){
    p$Feeding.Event.Link.Gap=Feeding.Event.Link.Gap
  }
  options(tmp.O)
  p
  
}

Get.Parameter.Names<-function(parameters){
  chambernames<-paste("Ch",1:length(parameters$Chamber.Sets),sep="")
  result<-c("BaselineWindowMin","FeedingThreshold","FeedingMinimum","TastingLow","TastingHigh",
            "FeedingMinEvents","TastingMinEvents","SamplesSec","ChamberSize",chambernames,"Link.Gap","PI.Multiplier")
  result
}

## This is used internally to return the group of parameters in vector form.
GetParameterVector<-function(parameters){  
  unlist(parameters)   
}
AreParametersEqual<-function(p1,p2){
  result<-TRUE
  if(p1$Baseline.Window.Minutes!=p2$Baseline.Window.Minutes) {
    result<-FALSE
  }
  
  if(p1$Feeding.Threshold!=p2$Feeding.Threshold) {
    result<-FALSE
  }
  if(p1$Feeding.Minimum!=p2$Feeding.Minimum) {
    result<-FALSE
  }
  if(p1$Tasting.Minimum!=p2$Tasting.Minimum) {
    result<-FALSE
  }
  if(p1$Tasting.Maximum!=p2$Tasting.Maximum) {
    result<-FALSE
  }
  if(p1$Feeding.Minevents!=p2$Feeding.Minevents){
    result<-FALSE
  }
  if(p1$Tasting.Minevents!=p2$Tasting.Minevents){
    result<-FALSE
  }
  if(p1$Chamber.Size!=p2$Chamber.Size) {
    result<-FALSE
  }
  if(p1$Feeding.Event.Link.Gap!=p2$Feeding.Event.Link.Gap) {
    result<-FALSE
  }
  if(p1$PI.Multiplier!=p2$PI.Multiplier) {
    result<-FALSE
  }
  result
}


### 3. DFMClass***************************************************************************

DFMClass<-function(id,parameters,range=c(0,0)){
  if (!is.numeric(id) || !all(is.finite(id)))
    stop("invalid arguments")
  ## Need to figure out which function to call based on
  ## DFM version and number of data files.
  
  ## First check if there are V3 files.
  tmp<-paste("DFM",id,"_.*[.]csv",sep="")
  files<-list.files(pattern=tmp)
  if(length(files)>0){
    DFMClassV3(id,parameters,range)
  }
  else {
    tmp<-paste("DFM_",id,".csv",sep="")
    files<-list.files(pattern=tmp)
    if(length(files)>0){
      DFMClassV2(id,parameters,range)
    }
    else {
      tmp<-paste("DFM_",id,"_.*[.]csv",sep="")
      files<-list.files(pattern=tmp)
      if(length(files)>0){
        DFMClassV2.LinkFiles(id,parameters,range)
      }
      else {
        print("DFM data file(s) not found.")  
      }
    }
  }
}



## These are functions that are generally available to the user ##
DFMClassV2<-function(id,parameters,range=c(0,0)) {
  if (!is.numeric(id) || !all(is.finite(id)))
    stop("invalid arguments")
  
  ## Check to determine whether the DFM object already exists
  st<-paste("DFM",id,sep="")
  found=0
  if(exists(st,where=1)) {
    data<-get(st)  
    found<-1
    if(AreParametersEqual(parameters,data$Parameters)==FALSE)
      data<-ChangeParameterObject(data,parameters)
  }
  
  ## If doesn't exist, get and create
  if(found==0) {
    
    file<-paste("DFM_",id,".csv",sep="")
    dfm<-read.csv(file,header=TRUE)  
    print(paste("Reading DFMV2 File:",file))
    ## Get Minutes from Sample column only if ElapsedTime is not
    ## there
    if('Seconds' %in% colnames(dfm)) {
      Minutes<-dfm$Seconds/60
      dfm<-data.frame(Minutes,dfm)
    } else if(('Date' %in% colnames(dfm))&&('Time' %in% colnames(dfm))&&('MSec' %in% colnames(dfm))){
      Seconds<-GetElapsedSeconds(dfm)
      Minutes<-Seconds/60.0
      dfm<-data.frame(Minutes,Seconds,dfm)
    } else {
      stop("Time information missing from DFM data.")
    }
    if(sum(range)!=0){
      dfm<- dfm[(dfm$Minutes>range[1]) & (dfm$Minutes<range[2]),]
    }
    data=list(ID=id,Parameters=parameters,RawData=dfm)
    class(data)="DFM"
    if(!is.na(FindDataBreaks(data,multiplier=4,returnvals=FALSE))){
      cat("Data lapses found. Use FindDataBreaks for details.")
      flush.console()      
    }
    data<-CalculateBaseline(data)  
    assign(st,data,pos=1)  
  }
  data 
}

DFMClassV3<-function(id,parameters,range=c(0,0)) {
  if (!is.numeric(id) || !all(is.finite(id)))
    stop("invalid arguments")
  
  ## Check to determine whether the DFM object already exists
  st<-paste("DFM",id,sep="")
  found=0
  if(exists(st,where=1)) {
    data<-get(st)  
    found<-1
    if(AreParametersEqual(parameters,data$Parameters)==FALSE)
      data<-ChangeParameterObject(data,parameters)
  }
  
  ## If doesn't exist, get and create
  if(found==0) {
    
    tmp<-paste("DFM",id,"_.*[.]csv",sep="")
    files<-list.files(pattern=tmp)
    files<-mixedsort(files)
    dfm<-read.csv(files[1],header=TRUE)  
    print(paste("Reading DFMV3 File:",files[1]))
    if(length(files)>1){
      for(i in 2:length(files)){ 
        print(paste("Reading DFMV3 File:",files[i]))
        tmp<-read.csv(files[i],header=TRUE)  
        dfm<-rbind(dfm,tmp)
      }
    }
    ## Get Minutes from Sample column only if ElapsedTime is not
    ## there
    if('Seconds' %in% colnames(dfm)) {
      Minutes<-dfm$Seconds/60
      dfm<-data.frame(Minutes,dfm)
    } else if(('Date' %in% colnames(dfm))&&('Time' %in% colnames(dfm))&&('MSec' %in% colnames(dfm))){
      Seconds<-GetElapsedSeconds(dfm)
      Minutes<-Seconds/60.0
      dfm<-data.frame(Minutes,Seconds,dfm)
    } else {
      stop("Time information missing from DFM data.")
    }
    if(sum(range)!=0){
      dfm<- dfm[(dfm$Minutes>range[1]) & (dfm$Minutes<range[2]),]
    }
    data=list(ID=id,Parameters=parameters,RawData=dfm)
    class(data)="DFM"
    if(!is.na(FindDataBreaks(data,multiplier=4,returnvals=FALSE))){
      cat("Data lapses found. Use FindDataBreaks for details.")
      flush.console()      
    }
    data<-CalculateBaseline(data)  
    assign(st,data,pos=1)  
  }
  data 
}

DFMClassV2.LinkFiles<-function(id,parameters,range=c(0,0)) {
  if (!is.numeric(id) || !all(is.finite(id)))
    stop("invalid arguments")
  
  ## Check to determine whether the DFM object already exists
  st<-paste("DFM",id,sep="")
  found=0
  if(exists(st,where=1)) {
    data<-get(st)  
    found<-1
    if(AreParametersEqual(parameters,data$Parameters)==FALSE)
      data<-ChangeParameterObject(data,parameters)
  }
  
  ## If doesn't exist, get and create
  if(found==0) {
    
    tmp<-paste("DFM_",id,"_.*[.]csv",sep="")
    files<-list.files(pattern=tmp)
    files<-mixedsort(files)
    dfm<-read.csv(files[1],header=TRUE)  
    print(paste("Reading DFM V2 File:",files[1]))
    if(length(files)>1){
      for(i in 2:length(files)){ 
        print(paste("Reading DFM V2 File:",files[i]))
        tmp<-read.csv(files[i],header=TRUE)  
        dfm<-rbind(dfm,tmp)
      }
    }
    ## Get Minutes from Sample column only if ElapsedTime is not
    ## there
    if('Seconds' %in% colnames(dfm)) {
      Minutes<-dfm$Seconds/60
      dfm<-data.frame(Minutes,dfm)
    } else if(('Date' %in% colnames(dfm))&&('Time' %in% colnames(dfm))&&('MSec' %in% colnames(dfm))){
      Seconds<-GetElapsedSeconds(dfm)
      Minutes<-Seconds/60.0
      dfm<-data.frame(Minutes,Seconds,dfm)
    } else {
      stop("Time information missing from DFM data.")
    }
    if(sum(range)!=0){
      dfm<- dfm[(dfm$Minutes>range[1]) & (dfm$Minutes<range[2]),]
    }
    data=list(ID=id,Parameters=parameters,RawData=dfm)
    class(data)="DFM"
    if(!is.na(FindDataBreaks(data,multiplier=4,returnvals=FALSE))){
      cat("Data lapses found. Use FindDataBreaks for details.")
      flush.console()      
    }
    data<-CalculateBaseline(data)  
    assign(st,data,pos=1)  
  }
  data 
}
## This function will look for consecutive entries in the
## RawData$Sec column whose difference is larger than it
## should be based on the Samples.Per.Sec parameter.
FindDataBreaks<-function(dfm,multiplier=4,returnvals=TRUE){
  Interval<-diff(dfm$RawData$Seconds)
  Interval<-c(0,Interval)
  thresh<-(1.0/dfm$Parameters$Samples.Per.Second)*multiplier
  Index<-1:length(Interval)
  Index<-Index[Interval>thresh]
  Interval<-Interval[Interval>thresh]
  if(returnvals==TRUE) {
    if(length(Interval)==0)
      c(NA)
    else
      cbind(Index,Interval,dfm$RawData[Index,])  
  }
  else {
    if(length(Interval)==0)
      c(NA)
    else
      c(1)
  } 
  
}

FindDataBreaksV3<-function(dfm,multiplier=4,returnvals=TRUE){
  diffs<-diff(dfm$RawData$Index)
  diffs<-c(0,diffs)
  thresh<-1
  Index<-1:length(diffs)
  Index<-Index[diffs>thresh]
  diffs<-diffs[diffs>thresh]
  if(returnvals==TRUE) {
    if(length(diffs)==0)
      c(NA)
    else
      cbind(Index,diffs,dfm$RawData[Index,])  
  }
  else {
    if(length(diffs)==0)
      c(NA)
    else
      c(1)
  } 
}

LastSampleData<-function(dfm){
  tmp<-BaselinedData(dfm)
  nr<-nrow(tmp)
  tmp[nr,]
}
FirstSampleData<-function(dfm){
  tmp<-BaselinedData(dfm)  
  tmp[1,]
}
BaselineData<-function(dfm,range=c(0,0)){
  tmp<-dfm$BaselineData
  if(sum(range)!=0) {
    tmp<- tmp[(tmp$Minutes>range[1]) & (tmp$Minutes<range[2]),]
  }    
  tmp  
}
RawData=function(dfm,range=c(0,0)) {
  tmp<-dfm$RawData
  if(sum(range)!=0) {
    tmp<- tmp[(tmp$Minutes>range[1]) & (tmp$Minutes<range[2]),]
  }    
  tmp
}
CleanDFM<-function(){
  tmp<-ls(pattern="DFM[0-9]",pos=1)
  rm(list=tmp,pos=1)
}
ChangeParameterObject<-function(dfm,newP) {
  p<-dfm$Parameters
  baseline.flag<-FALSE
  threshold.flag<-FALSE
  eventpi.flag<-FALSE
  tmp.O<-options()
  options(warn=-1)
  dfm$Parameters<-newP
  ## Change only those that are listed
  if(p$Baseline.Window.Minutes!=newP$Baseline.Window.Minutes) {    
    baseline.flag<-TRUE
  }
  if(p$Feeding.Threshold!=newP$Feeding.Threshold) {    
    threshold.flag<-TRUE
  }
  if(p$Feeding.Minimum!=newP$Feeding.Minimum) {    
    threshold.flag<-TRUE
  }
  if(p$Tasting.Minimum!=newP$Tasting.Minimum) {    
    threshold.flag<-TRUE
  }
  if(p$Tasting.Maximum!=newP$Tasting.Maximum) {    
    threshold.flag<-TRUE
  }
  if(p$Feeding.Minevents!=newP$Feeding.Minevents){
    eventpi.flag<-TRUE
  }
  if(p$Tasting.Minevents!=newP$Tasting.Minevents){
    eventpi.flag<-TRUE
  }
  if(p$Samples.Per.Second!=newP$Samples.Per.Second){
    adaptive.baseline.flag<-TRUE
  }
  if(p$Chamber.Size !=newP$Chamber.Size){
    baseline.flag<-TRUE
  }
  if(p$Feeding.Event.Link.Gap != newP$Feeding.Event.Link.Gap){
    threshold.flag<-TRUE
  }
  if(sum(c(p$Chamber.Sets)!=c(newP$Chamber.Sets))!=0){
    baseline.flag<-TRUE
  }
  
  if(p$PI.Multiplier!=newP$PI.Multiplier){
    eventpi.flag<-TRUE
  }
  
  ## Now update the stats needed
  if(baseline.flag==TRUE) {
    dfm<-CalculateBaseline(dfm)
  }
  else if(threshold.flag==TRUE) {
    dfm<-SetThreshold(dfm,getStandard=FALSE)
  }
  else if(eventpi.flag==TRUE) {
    dfm<-Set.Feeding.Data(dfm)
    dfm<-Set.Tasting.Data(dfm)
    dfm<-Set.Durations.And.Intervals(dfm)
    dfm<-Set.Tasting.Durations.And.Intervals(dfm)
  }
  options(tmp.O)
  UpdateHiddenDFMObject(dfm)
  dfm
}

GetDFM<-function(id){
  if (!is.numeric(id) || !all(is.finite(id)))
    stop("invalid arguments")
  
  ## Check to determine whether the DFM object already exists
  st<-paste("DFM",id,sep="")
  data<-NA
  if(exists(st,where=1)) {
    data<-get(st)  
  }
  data
}

#**************************************************************************************
#*
#*
#*
#*
#*
#  MiscFunctionsR



## Give filename to this function and it will do analyses of 
## variance with and without DFM as a factor.  It will also plot
## the data for each well for each treatment in each DFM to look for
## consistency.
QuickAOV.FeedingSummary<-function(datafile="FeedingSummary.csv"){
  data<-read.csv(datafile)
  data$Licks<-data$Licks^0.25
  data$DFM<-factor(data$DFM)
  
  aov.1<-aov(Licks~Treatment,data=data)
  print(summary(aov.1))
  cat("\n********\n")
  aov.2<-aov(Licks~DFM+Treatment,data=data)
  print(summary(aov.2))
  
  g1<-ggplot(data,aes(x=DFM,y=Licks)) + geom_boxplot() + geom_jitter(aes(color=Treatment))
  show(g1)
  
  
  g2<-ggplot(data,aes(x=DFM,y=Licks, fill=Treatment)) + geom_boxplot(aes(fill=Treatment))
  show(g2)
}
PlotLicksandLight.Well<-function(dfm,well,range=c(0,0),TransformLicks=TRUE){
  tmp<-FeedingData.Well.Licks(dfm,well)
  SumLicks<-cumsum(tmp)
  if(TransformLicks==TRUE)
    SumLicks<-SumLicks^0.25
  Lights<-GetLightsInfo(dfm)
  Light<-Lights[,paste("W",well,sep="")]
  row<-1
  col<-1
  Row<-rep(row,length(tmp))
  Col<-rep(col,length(tmp))
  Well<-rep(1,length(tmp))
  Size <- Light*2
  Size[Size==0]<-1
  results<-data.frame(dfm$LickData$Minutes,SumLicks,Light,Size,Row,Col,Well)
  names(results)<-c("Minutes","SumLicks","Light","Size","Row","Col","Well")
  if(sum(range)!=0) {
    results<- results[(results$Minutes>range[1]) & (results$Minutes<range[2]),]
  }   
  if(TransformLicks==TRUE)
    ylabel="Transformed Cumulative Licks"
  else
    ylabel="Cumulative Licks"
  gp<-ggplot(results,aes(Minutes,SumLicks,color=Light)) + geom_point(size=results$Size) +
    ggtitle(paste("DFM ",dfm$ID, "; Well ",well, sep=""))+ ylab(ylabel)+ labs(color="Light")
  
  gp
}


GetLightsInfo<-function(dfm, range=c(0,0)){
  rd<-RawData(dfm,range)
  Column1<-rd$OptoCol1
  Column2<-rd$OptoCol2
  data<-data.frame(Column1,Column2)
  data2<-apply(data,1,LightStatus)
  data2<-t(data2)
  final.data<-data.frame(rd$Minutes,data2,data)
  names(final.data)<-c("Minutes",paste("W",1:12,sep=""),"OptoCol1","OptoCol2")
  final.data
}


LightStatus<-function(cols){
  col1<-cols[1]
  col2<-cols[2]
  lights<-rep(FALSE,12)
  for(i in 0:5){
    tmp<-bitwShiftL(1,i)
    if(bitwAnd(tmp,col1)>0){
      lights[i*2+1]<-TRUE
    }
    if(bitwAnd(tmp,col2)>0){
      lights[i*2+2]<-TRUE
    }
  }
  lights
}


ScrollEventPlots<-function(dfm,wellnum){
  events<-unique(dfm$BaselineData$Sample)
  currentevent<-1
  keepgoing<-TRUE
  while (keepgoing==TRUE){
    cat("Return = Forward; b = back; q = quit")
    PlotSingleSampleEvent(dfm,wellnum,events[currentevent])
    tmp<-readline()
    if(tmp[1]==""){
      currentevent <- currentevent+1
      if(currentevent>length(events))
        currentevent<-1
    }
    else if(tmp[1]=="b"){
      currentevent<-currentevent-1
      if(currentevent<1)
        currentevent=length(events)
    }
    
    else if(tmp[1]=="q"){
      keepgoing<-FALSE
    }
  }
}


PlotSingleSampleEvent<-function(dfm,wellnum,eventnum){
  tmp<-subset(dfm$BaselineData,dfm$BaselineData$Sample==eventnum)
  w<-paste("W",wellnum,sep="")
  gp<-ggplot(tmp,aes(Minutes,tmp[,(6+wellnum)])) + geom_line(color="red",size=1.2) + facet_wrap(tmp$Sample) +geom_point(color="blue",size=4) +
    ggtitle(paste("DFM: ",dfm$ID,"   Well: ",w,"   Event:",eventnum)) + ylab("Signal") + labs(color="Chamber")
  show(gp)
}



GetCircMin<-function(a){
  
}

# THis accepts the startTime as MM/DD/YYYY HH:MM:SS in military time
# The data frame argument must have a column entitled "Minutes", which
# is assumed to be the elapsed time of the experiment.
AddRealTimeColumns<-function(data,startTimeString){
  dtm<-strptime(startTimeString, format = "%m/%d/%Y %H:%M:%S")
  secs<-data$Minutes*60
  realTimes<-as.POSIXlt(dtm+secs)
  tmp<-names(data)
  cm<-realTimes$hour*60+realTimes$min+realTimes$sec/60
  tmp<-c("RealTime","CircMinutes",tmp)
  r<-data.frame(realTimes,cm,data)
  names(r)<-tmp
  r
}


## Beta Code for examining summary feeding parameters by a fixed number of
## events rather than a fixed time interval.
GetEventRangeInMinutes.Well.OneWell<-function(dfm,well,events=c(0,0)){
  tmp<-FeedingData.Events(dfm)
  cname=paste("W",well,sep="")
  data<-tmp[,cname]  
  min<-tmp$Minutes
  data<-cumsum(data>0)
  if(sum(events)!=0) {
    start<-events[1]
    end<-events[2]
    #print(data[data>=start & data<=end])
    min2<-min[data>=start & data<=end]
  }
  else {
    min2<-min
  }
  if(length(min2)==0){
    result<-c(min[length(min)],min[length(min)])  
  }
  else {
    result<-c(min2[1],min2[length(min2)])
  }
  
  result
}
GetEventRangeInMinutes.Well.TwoWell<-function(dfm,chamber,events=c(0,0)){
  tmp<-FeedingData.Events(dfm)
  if(chamber==1){
    d1<-tmp[,"W1"]
    d2<-tmp[,"W2"]
  }else if(chamber==2){
    d1<-tmp[,"W3"]
    d2<-tmp[,"W4"]
  }else if(chamber==3){
    d1<-tmp[,"W5"]
    d2<-tmp[,"W6"]
  }else if(chamber==4){
    d1<-tmp[,"W7"]
    d2<-tmp[,"W8"]
  }else if(chamber==5){
    d1<-tmp[,"W9"]
    d2<-tmp[,"W10"]
  }else if(chamber==6){
    d1<-tmp[,"W11"]
    d2<-tmp[,"W12"]
  }
  data<-d1+d2
  min<-tmp$Minutes
  data<-cumsum(data>0)
  if(sum(events)!=0) {
    start<-events[1]
    end<-events[2]
    #print(data[data>=start & data<=end])
    min2<-min[data>=start & data<=end]
  }
  else {
    min2<-min
  }
  if(length(min2)==0){
    result<-c(min[length(min)],min[length(min)])  
  }
  else {
    result<-c(min2[1],min2[length(min2)])
  }
  
  result
}
Feeding.Summary.TwoWell.ByEvents<-function(dfm,events=c(0,0),TransformLicks=TRUE){
  if(dfm$Parameters$Chamber.Size!=2)
    stop("This function is for two-chamber DFM only")
  
  for(i in 1:nrow(dfm$Parameters$Chamber.Sets)) {
    range<-GetEventRangeInMinutes.Well.TwoWell(dfm,i,events)
    # Need a small correction here to avoid missing the first event
    range[1]<-range[1]-0.0034
    
    
    if(dfm$Parameters$PI.Multiplier==1){
      wellA<-dfm$Parameters$Chamber.Sets[i,1]
      wellB<-dfm$Parameters$Chamber.Sets[i,2] 
    }
    else {
      wellB<-dfm$Parameters$Chamber.Sets[i,1]
      wellA<-dfm$Parameters$Chamber.Sets[i,2] 
    }
    
    interval.a<-Feeding.IntervalSummary.Well(dfm,wellA,range)
    intensity.a<-Feeding.IntensitySummary.Well(dfm,wellA,range)
    dur.a<-Feeding.DurationSummary.Well(dfm,wellA,range)  
    FLicks.a<-Feeding.TotalLicks.Well(dfm,wellA,range)
    FEvents.a<-Feeding.TotalEvents.Well(dfm,wellA,range)    
    
    interval.b<-Feeding.IntervalSummary.Well(dfm,wellB,range)
    intensity.b<-Feeding.IntensitySummary.Well(dfm,wellB,range)
    dur.b<-Feeding.DurationSummary.Well(dfm,wellB,range)  
    FLicks.b<-Feeding.TotalLicks.Well(dfm,wellB,range)
    FEvents.b<-Feeding.TotalEvents.Well(dfm,wellB,range) 
    
    FPIs<-c((FLicks.a-FLicks.b)/(FLicks.a+FLicks.b),(FEvents.a-FEvents.b)/(FEvents.a+FEvents.b))
    if(i==1){
      result<-data.frame(matrix(c(dfm$ID,i,FPIs,FLicks.a,FLicks.b,FEvents.a,FEvents.b,unlist(dur.a),unlist(dur.b),unlist(interval.a),
                                  unlist(interval.b),unlist(intensity.a),unlist(intensity.b),range[1],range[2]),nrow=1))    
      
    }
    else {
      tmp<-data.frame(matrix(c(dfm$ID,i,FPIs,FLicks.a,FLicks.b,FEvents.a,FEvents.b,unlist(dur.a),unlist(dur.b),unlist(interval.a),
                               unlist(interval.b),unlist(intensity.a),unlist(intensity.b),range[1],range[2]),nrow=1))
      result<-rbind(result,tmp)      
    }
  }
  names(result)<-c("DFM","Chamber","PI","EventPI","LicksA","LicksB","EventsA","EventsB","MeanDurationA","MedDurationA",
                   "MeanDurationB","MedDurationB","MeanTimeBtwA","MedTimeBtwA",
                   "MeanTimeBtwB","MedTimeBtwB","MeanIntA","MedianIntA","MinIntA","MaxIntA",
                   "MeanIntB","MedianIntB","MinIntB","MaxIntB","StartMin","EndMin")
  if(TransformLicks==TRUE){
    result$LicksA<-result$LicksA^0.25
    result$LicksB<-result$LicksB^0.25
  }
  result    
  
}
Feeding.Summary.OneWell.ByEvents<-function(dfm,events=c(0,0),TransformLicks=TRUE){
  if(dfm$Parameters$Chamber.Size!=1)
    stop("This function is for single chambers only")
  
  for(i in 1:12 ){
    range<-GetEventRangeInMinutes.Well.OneWell(dfm,i,events)
    # Need a small correction here to avoid missing the first event
    range[1]<-range[1]-0.0034
    interval<-Feeding.IntervalSummary.Well(dfm,i,range)
    intensity<-Feeding.IntensitySummary.Well(dfm,i,range)
    dur<-Feeding.DurationSummary.Well(dfm,i,range)  
    FLicks<-Feeding.TotalLicks.Well(dfm,i,range)
    FEvents<-Feeding.TotalEvents.Well(dfm,i,range)    
    if(i==1)
      result<-data.frame(matrix(c(dfm$ID,i,FLicks,FEvents,unlist(dur),unlist(interval),unlist(intensity),range[1],range[2]),nrow=1))  
    else {
      tmp<-data.frame(matrix(c(dfm$ID,i,FLicks,FEvents,unlist(dur),unlist(interval),unlist(intensity),range[1],range[2]),nrow=1))  
      result<-rbind(result,tmp)
    }      
  }
  names(result)<-c("DFM","Chamber","Licks","Events","MeanDuration","MedDuration",
                   "MeanTimeBtw","MedTimeBtw","MeanInt","MedianInt","MinInt","MaxInt","StartMin","EndMin")
  if(TransformLicks==TRUE)
    result$Licks<-result$Licks^0.25
  result    
}
Feeding.Summary.DFM.ByEvents<-function(dfm,events=c(0,0),TransformLicks=TRUE){
  if(dfm$Parameters$Chamber.Size==1)
    Feeding.Summary.OneWell.ByEvents(dfm,events,TransformLicks)
  else if(dfm$Parameters$Chamber.Size==2)
    Feeding.Summary.TwoWell.ByEvents(dfm,events,TransformLicks)
  else
    stop("Feeding Summary not implemented for this DFM type.")    
}

Feeding.Summary.DFM.BySingleEvents<-function(dfm,events=c(1,2),TransformLicks=TRUE){
  events<-events[1]:events[2]
  for(i in 1:(length(events))){
    new.events<-c(events[i],events[i]+0.5)
    if(dfm$Parameters$Chamber.Size==1)
      tmp<-Feeding.Summary.OneWell.ByEvents(dfm,new.events,TransformLicks)
    else if(dfm$Parameters$Chamber.Size==2)
      tmp<-Feeding.Summary.TwoWell.ByEvents(dfm,new.events,TransformLicks)
    else
      stop("Feeding Summary not implemented for this DFM type.")    
    EventNum<-rep(events[i],nrow(tmp))
    tmp<-data.frame(EventNum,tmp)
    if(i==1){
      results<-tmp
    }
    else {
      results<-rbind(results,tmp)
    }
  }
  results
}

Feeding.Summary.Monitors.BySingleEvents<-function(monitors,parameters,events=c(1,2),expDesign=NA,SaveToFile=TRUE,TransformLicks=TRUE,filename="FeedingSummaryBySingleEvents"){
  individ.params<-FALSE
  ## Check to determine whether parameters is a signle parameter object
  ## or a list of them.  If it is a single one, then we use the same one for all
  ## if it is a list, then we use a different one for each.
  if(is.list(parameters[[1]])==TRUE){
    if(length(parameters)!=length(monitors))
      stop("If individuals parameter objects are specified, there must be one for each DFM.")
    individ.params<-TRUE
  }
  
  for(j in 1:length(monitors)){
    monitor<-monitors[j]
    if(individ.params==TRUE)
      p<-parameters[[j]]
    else
      p<-parameters
    dfm<-DFMClass(monitor,p)  
    parameter.vector<-matrix(GetParameterVector(p),nrow=1)
    pnames<-Get.Parameter.Names(p)
    tmp<-Feeding.Summary.DFM.BySingleEvents(dfm,events,TransformLicks)      
    tmp2<-data.frame(tmp,parameter.vector)
    names(tmp2)<-c(names(tmp),pnames)
    if(j==1){
      results<-tmp2
    }
    else {
      results<-rbind(results,tmp2)  
    }
    
  }  
  if(is.data.frame(expDesign)) {
    results<-AppendTreatmentonResultsFrame(results,expDesign)
    trt.summary<-suppressWarnings(AggregateTreatments(results))
    if(SaveToFile==TRUE){
      filename2<-paste(filename,"_Stats.csv",sep="")
      write.csv(trt.summary,file=filename2,row.names=FALSE)
      filename2<-paste(filename,"_Data.csv",sep="")
      write.csv(results,file=filename2,row.names=FALSE)
    }
  }
  else if(SaveToFile==TRUE){
    filename<-paste(filename,".csv",sep="")
    write.csv(results,file=filename,row.names=FALSE)
  }
  
  if(is.data.frame(expDesign)) {
    return(list(Results=results,Stats=trt.summary))
  }
  else {
    return(list(Results=results))
  }
}


BinnedFeeding.Summary.DFM.ByEvents<-function(dfm,binsize.Enum,TransformLicks=TRUE){
  tmp<-Feeding.TotalEvents(dfm)
  y<-seq(0,max(tmp),by=binsize.Enum)
  
  tmpMatrix<-cbind(y[-length(y)],y[-1])
  tmpMatrix[,1]<-tmpMatrix[,1]+1
  intervals<-cut(y+0.000001,y,include.lowest=TRUE,dig.lab=8)
  intervals<-intervals[-length(intervals)]
  result<-Feeding.Summary.DFM.ByEvents(dfm,events=tmpMatrix[1,],TransformLicks)
  Interval<-rep(intervals[1],nrow(result))
  Minutes<-rep(mean(tmpMatrix[1,]),nrow(result))
  result<-data.frame(Interval,Minutes,result)
  for(i in 2:nrow(tmpMatrix)){
    tmp<-Feeding.Summary.DFM.ByEvents(dfm,events=tmpMatrix[i,],TransformLicks)
    Interval<-rep(intervals[i],nrow(tmp))
    Minutes<-rep(mean(tmpMatrix[i,]),nrow(tmp))
    tmp<-data.frame(Interval,Minutes,tmp)
    result<-rbind(result,tmp)
  }
  result
}
BinnedFeeding.Summary.Monitors.ByEvents<-function(monitors,parameters,binsize.Enum=5,expDesign=NA,SaveToFile=TRUE,TransformLicks=TRUE,filename="BinnedSummaryByEvents"){
  individ.params<-FALSE
  ## Check to determine whether parameters is a signle parameter object
  ## or a list of them.  If it is a single one, then we use the same one for all
  ## if it is a list, then we use a different one for each.
  if(is.list(parameters[[1]])==TRUE){
    if(length(parameters)!=length(monitors))
      stop("If individuals parameter objects are specified, there must be one for each DFM.")
    individ.params<-TRUE
  }
  
  for(j in 1:length(monitors)){
    monitor<-monitors[j]
    if(individ.params==TRUE)
      p<-parameters[[j]]
    else
      p<-parameters
    dfm<-DFMClass(monitor,p)  
    parameter.vector<-matrix(GetParameterVector(p),nrow=1)
    pnames<-Get.Parameter.Names(p)
    tmp<-BinnedFeeding.Summary.DFM.ByEvents(dfm,binsize.Enum,TransformLicks)      
    tmp2<-data.frame(tmp,parameter.vector)
    names(tmp2)<-c(names(tmp),pnames)
    if(j==1){
      results<-tmp2
    }
    else {
      results<-rbind(results,tmp2)  
    }
    
  }  
  if(is.data.frame(expDesign)) {
    results<-AppendTreatmentonResultsFrame(results,expDesign)
    trt.summary<-suppressWarnings(AggregateTreatmentsBinnedData(results))
    if(SaveToFile==TRUE){
      filename2<-paste(filename,"_Stats.csv",sep="")
      write.csv(trt.summary,file=filename2,row.names=FALSE)
      filename2<-paste(filename,"_Data.csv",sep="")
      write.csv(results,file=filename2,row.names=FALSE)
    }
  }
  else if(SaveToFile==TRUE){
    filename<-paste(filename,".csv",sep="")
    write.csv(results,file=filename,row.names=FALSE)
  }
  
  if(is.data.frame(expDesign)) {
    return(list(Results=results,Stats=trt.summary))
  }
  else {
    return(list(Results=results))
  }
}


## ****************************************************************
# FLICUserFunctions

##### QC Functions ######
OutputBaselinedData.DFM<-function(dfm, range=c(0,0),filename="Baselined"){
  
  parameter.vector<-GetDFMParameterVector(dfm)
  pnames<-Get.Parameter.Names(dfm$Parameters)
  
  tmp.all<-BaselineData(dfm,range)
  
  tmp<-data.frame(cbind(rep(dfm$ID,nrow(tmp.all)),tmp.all))
  tmp2<-matrix(rep(parameter.vector,nrow(tmp)),ncol=length(parameter.vector),byrow=TRUE)
  
  tmp3<-cbind(tmp,tmp2)
  final.names<-c("DFM",names(tmp.all),pnames)
  names(tmp3)<-final.names
  
  filename<-paste(filename,"_DFM",dfm$ID,".csv",sep="") 
  write.csv(tmp3,file=filename,row.names=FALSE)  
}
RawDataPlot.DFM<-function(dfm,range=c(0,0),OutputPNGFile=FALSE) {
  ##windows(record=FALSE,width=8,height=12) # opens a window and starts recording
  tmp<-RawData(dfm,range)
  tmp<-tmp[,c(1,7,8,9,10,11,12,13,14,15,16,17,18)]
  thresh.high<-dfm$Parameters$Feeding.Threshold
  thresh.low<-dfm$Parameters$Feeding.Minimum
  tmp.m<-melt(tmp, id="Minutes")
  gp<-ggplot(data = tmp.m, aes(x = Minutes, y = value)) +
    geom_line() + facet_grid(variable ~ .) + ggtitle(paste("DFM:",dfm$ID))
  fn<-paste("RawDataPlot_DFM",dfm$ID,".png",sep="")
  if(OutputPNGFile==TRUE)
    ggsave(filename = fn)
  else
    show(gp)
}
BaselinedDataPlot.DFM<-function(dfm,range=c(0,0),OutputPNGFile=FALSE, IncludeThresholds=FALSE) {
  ##windows(record=FALSE,width=8,height=12) # opens a window and starts recording
  tmp<-BaselineData(dfm,range)
  tmp<-tmp[,c(1,7,8,9,10,11,12,13,14,15,16,17,18)]
  thresh.high<-dfm$Parameters$Feeding.Threshold
  thresh.low<-dfm$Parameters$Feeding.Minimum
  tmp.m<-melt(tmp, id="Minutes")
  gp<-ggplot(data = tmp.m, aes(x = Minutes, y = value)) +
    geom_line() + facet_grid(variable ~ .) + ggtitle(paste("DFM:",dfm$ID))
  if(IncludeThresholds==TRUE)
    gp<-gp+
    geom_hline(yintercept=thresh.high,linetype="dashed", color = "red", size=0.2) + 
    geom_hline(yintercept=thresh.low,linetype="dashed", color = "red", size=0.2) 
  fn<-paste("BaselinedDataPlot_DFM",dfm$ID,".png",sep="")
  if(OutputPNGFile==TRUE)
    ggsave(filename = fn)
  else
    show(gp)
}
## Use this function to identify large well- or DFM-based effects.
PlotTotalLicks.Monitors<-function(monitors, p, range=c(0,0),TransformLicks=TRUE){
  tmp2<-Feeding.Summary.Monitors(monitors,p,range=range,TransformLicks)$Results
  
  ylabel="Licks"
  if(TransformLicks==TRUE) {
    ylabel="Transformed Licks"
  }
  
  tmp2$DFM<-factor(tmp2$DFM)
  
  if("LicksA" %in% names(tmp2)){
    tmp2<-tmp2[,c("DFM","LicksA","LicksB")]
    tmp2<-melt(tmp2,id.vars="DFM",value.name="Licks", variable.name="Well")
    gp<-ggplot(tmp2,aes(x=DFM,y=Licks,fill=Well)) + geom_dotplot(binaxis='y',stackdir='center',stackratio=1.5, dotsize=0.7) +
      ylab(ylabel)
  }
  else {
    gp<-ggplot(tmp2,aes(x=DFM,y=Licks,fill=DFM)) + geom_dotplot(binaxis='y',stackdir='center',stackratio=1.5, dotsize=0.7) +
      ylab(ylabel)
  }
  gp
}

##### Experiment Analysis Functions ######
Feeding.Summary.Monitors<-function(monitors,parameters,expDesign=NA,range=c(0,0),SaveToFile=TRUE,TransformLicks=TRUE,filename="FeedingSummary"){
  individ.params<-FALSE
  ## Check to determine whether parameters is a signle parameter object
  ## or a list of them.  If it is a single one, then we use the same one for all
  ## if it is a list, then we use a different one for each.
  if(is.list(parameters[[1]])==TRUE){
    if(length(parameters)!=length(monitors))
      stop("If individuals parameter objects are specified, there must be one for each DFM.")
    individ.params<-TRUE
  }
  
  for(j in 1:length(monitors)){
    ##print(paste("Summarizing feeding data for DFM ",monitors[j],".",sep=""))
    ##flush.console()
    monitor<-monitors[j]
    if(individ.params==TRUE)
      p<-parameters[[j]]
    else
      p<-parameters
    dfm<-DFMClass(monitor,p)  
    parameter.vector<-matrix(GetParameterVector(p),nrow=1)
    pnames<-Get.Parameter.Names(p)
    tmp<-Feeding.Summary.DFM(dfm,range,TransformLicks)      
    tmp2<-data.frame(tmp,parameter.vector)
    names(tmp2)<-c(names(tmp),pnames)
    if(j==1){
      results<-tmp2
    }
    else {
      results<-rbind(results,tmp2)  
    }
    
  }  
  if(is.data.frame(expDesign)) {
    results<-AppendTreatmentonResultsFrame(results,expDesign)
    trt.summary<-suppressWarnings(AggregateTreatments(results))
    if(SaveToFile==TRUE){
      filename2<-paste(filename,"_Stats.csv",sep="")
      write.csv(trt.summary,file=filename2,row.names=FALSE)
      filename2<-paste(filename,"_Data.csv",sep="")
      write.csv(results,file=filename2,row.names=FALSE)
    }
  }
  else if(SaveToFile==TRUE){
    filename<-paste(filename,".csv",sep="")
    write.csv(results,file=filename,row.names=FALSE)
  }
  
  if(is.data.frame(expDesign)) {
    return(list(Results=results,Stats=trt.summary))
  }
  else {
    return(list(Results=results))
  }
}
BinnedFeeding.Summary.Monitors<-function(monitors,parameters,binsize.min=30,expDesign=NA,range=c(0,0),SaveToFile=TRUE,TransformLicks=TRUE,filename="BinnedSummary"){
  individ.params<-FALSE
  ## Check to determine whether parameters is a signle parameter object
  ## or a list of them.  If it is a single one, then we use the same one for all
  ## if it is a list, then we use a different one for each.
  if(is.list(parameters[[1]])==TRUE){
    if(length(parameters)!=length(monitors))
      stop("If individuals parameter objects are specified, there must be one for each DFM.")
    individ.params<-TRUE
  }
  
  ## If there is no range specified, choose the last time for all monitors.
  if(sum(range)==0) {
    EndTimeMin<-0
    for(j in 1:length(monitors)){
      monitor<-monitors[j]
      if(individ.params==TRUE)
        p<-parameters[[j]]
      else
        p<-parameters
      dfm<-DFMClass(monitor,p)   
      tmp<-LastSampleData(dfm)$Minutes
      if(tmp>EndTimeMin)
        EndTimeMin<-tmp
    }
    range[2]<-EndTimeMin
  }
  
  for(j in 1:length(monitors)){
    monitor<-monitors[j]
    if(individ.params==TRUE)
      p<-parameters[[j]]
    else
      p<-parameters
    dfm<-DFMClass(monitor,p)  
    parameter.vector<-matrix(GetParameterVector(p),nrow=1)
    pnames<-Get.Parameter.Names(p)
    tmp<-BinnedFeeding.Summary.DFM(dfm,binsize.min,range,TransformLicks)      
    tmp2<-data.frame(tmp,parameter.vector)
    names(tmp2)<-c(names(tmp),pnames)
    if(j==1){
      results<-tmp2
    }
    else {
      results<-rbind(results,tmp2)  
    }
    
  }  
  if(is.data.frame(expDesign)) {
    results<-AppendTreatmentonResultsFrame(results,expDesign)
    trt.summary<-suppressWarnings(AggregateTreatmentsBinnedData(results))
    if(SaveToFile==TRUE){
      filename2<-paste(filename,"_Stats.csv",sep="")
      write.csv(trt.summary,file=filename2,row.names=FALSE)
      filename2<-paste(filename,"_Data.csv",sep="")
      write.csv(results,file=filename2,row.names=FALSE)
    }
  }
  else if(SaveToFile==TRUE){
    filename<-paste(filename,".csv",sep="")
    write.csv(results,file=filename,row.names=FALSE)
  }
  
  if(is.data.frame(expDesign)) {
    return(list(Results=results,Stats=trt.summary))
  }
  else {
    return(list(Results=results))
  }
}

##### Simple Experiment Plotting Functions ######
## These functions take the output of the experimental analysis functions and rapidly plot different measures.
## Current Type options are (case sensitive): Licks, Events, Durations, MinInt, TimeBtw, PI, EventPI (latter 2 for choice chambers only).
BinnedDataPlot<-function(binnedDataResult,Type="Licks",SaveToFile=FALSE){
  if("LicksA" %in% names(binnedDataResult$Results)){
    BinnedPlot.TwoWell.Trt(binnedDataResult,Type,SaveToFile)
  }
  else {
    BinnedPlot.OneWell.Trt(binnedDataResult,Type,SaveToFile)
  }
}

## Current Type options are (case sensitive): Licks, Events, Durations, MinInt, TimeBtw, PI, EventPI (latter 2 for choice chambers only).
DataPlot<-function(summaryResults,Type="Licks",SaveToFile=FALSE){
  if("LicksA" %in% names(summaryResults$Results)){
    SimpleDataPlot.TwoWell(summaryResults,Type,SaveToFile)
  }
  else {
    SimpleDataPlot.OneWell(summaryResults,Type,SaveToFile)
  }
}

##### Advanced Cumulative Division Box Plots ######
## Divisions will create a separate graph for the cumulative PI (starting at range[0]) up
## to each of the division points in the experiment.  This is distinct from the time-dependent
## PI because it will always start from the first point in the data set (satisfying the first
## range parameter).

## These plots allow the user to specify a range.

## Although it is not optimal, the two-well aspect of these functions just adds or averages the values of
## both wells.  In the future this should probably present values for each well, adjusted for the PI multiplier.

## These plots are slower than the simple plots because analysis calculations are required.
## Division plots are not optimized for choice chambers because it combines well data at the moment (4/23/2020).
## Current Type options are (case sensitive): Licks, Events, Durations, MinInt, TimeBtw, PI, EventPI (latter 2 for choice chambers only).
DivisionPlots<-function(monitors,parameters,expDesign,range=c(0,0),divisions=1,Type="Licks",SaveToFile=FALSE,TransformLicks=TRUE){
  if(is.list(parameters[[1]])){
    cs <- parameters[[1]]$Chamber.Size
  }
  else {
    cs<-parameters$Chamber.Size
  }
  if(cs==1)
    DivisionPlots.OneWell(monitors,parameters,expDesign,range,divisions,Type,SaveToFile,TransformLicks)
  else if(cs==2)
    DivisionPlots.TwoWell(monitors,parameters,expDesign,range,divisions,Type,SaveToFile,TransformLicks)
  else
    stop("Feeding lick plots not implemented for this DFM type.")    
}

##### Cumulative plots for choice chambers only (Ugly) #####
CumulativePIPlots<-function(monitors,parameters,expDesign,range=c(0,0),SaveToFile=FALSE){
  individ.params<-FALSE
  if(is.list(parameters[[1]])==TRUE){
    if(length(parameters)!=length(monitors))
      stop("If individuals parameter objects are specified, there must be one for each DFM.")
    individ.params<-TRUE
  }
  
  if(SaveToFile==TRUE){
    filename<-paste("CumulativePI.pdf",sep="")
    pdf(file=filename)
  }
  
  ## Basically plot the culuative PI plot for each chamber for each treatment
  ## Each Treatment will get a different color.
  for(i in 1:length(monitors)){
    if(individ.params==TRUE)
      p<-parameters[[i]]
    else
      p<-parameters
    dfm<-DFMClass(monitors[i],p)
    tmp<-CumulativePI.DFM(dfm,ShowPlots = FALSE)
    DFM<-rep(dfm$ID,nrow(tmp))
    tmp<-data.frame(tmp,DFM)
    if(!exists("results",inherits = FALSE))
      results<-tmp
    else
      results<-rbind(results,tmp)
    
  }
  Treatment<-apply(results,1,GetTreatmentForRow,expDesign)
  
  results<-data.frame(results,Treatment)
  results<-subset(results,Treatment!="None")
  
  ## Note that we can just average the curves here because different monitors will have slightly different
  ## values for the minutes column.  They were not all collected at exactly the same time.  So aggregating
  ## won't work.
  p<-ggplot(results, aes(Minutes, PI, group=interaction(DFM,Chamber))) + geom_line(aes(color=Treatment),size=1) +
    geom_smooth(aes(x=Minutes,y = PI,group = Treatment,color=Treatment),size=3)
  show(p)
  
  if(SaveToFile==TRUE){ 
    graphics.off()
  }
  
}
CumulativeEventPIPlots<-function(monitors,parameters,expDesign,events.limit=NA,by.bout=FALSE,SaveToFile=FALSE){
  individ.params<-FALSE
  if(is.list(parameters[[1]])==TRUE){
    if(length(parameters)!=length(monitors))
      stop("If individuals parameter objects are specified, there must be one for each DFM.")
    individ.params<-TRUE
  }
  
  if(SaveToFile==TRUE){
    filename<-paste("CumulativeEventPI.pdf",sep="")
    pdf(file=filename)
  }
  
  ## Basically plot the culuative PI plot for each chamber for each treatment
  ## Each Treatment will get a different color.
  for(i in 1:length(monitors)){
    if(individ.params==TRUE)
      p<-parameters[[i]]
    else
      p<-parameters
    dfm<-DFMClass(monitors[i],p)
    tmp<-CumulativeEventPI.DFM(dfm,events.limit,ShowPlots = FALSE)
    DFM<-rep(dfm$ID,nrow(tmp))
    tmp<-data.frame(tmp,DFM)
    if(!exists("results",inherits = FALSE))
      results<-tmp
    else
      results<-rbind(results,tmp)
    
  }
  Treatment<-apply(results,1,GetTreatmentForRow,expDesign)
  
  results<-data.frame(results,Treatment)
  results<-subset(results,Treatment!="None")
  
  ## Note that we have to smooth the curves here because different monitors will have slightly different
  ## values for the minutes column.  They were not all collected at exactly the same time.  So aggregating
  ## won't work.
  if(by.bout==TRUE){
    p<-ggplot(results, aes(EventNum, PI, group=interaction(DFM,Chamber))) + geom_line(aes(color=Treatment),size=1) +
      geom_smooth(aes(x=EventNum,y = PI,group = Treatment,color=Treatment),size=3) + xlab("Event Number")+ ylab("PI (Events)")
  }
  else {
    p<-ggplot(results, aes(Minutes, PI, group=interaction(DFM,Chamber))) + geom_line(aes(color=Treatment),size=1) +
      geom_smooth(aes(x=Minutes,y = PI,group = Treatment,color=Treatment),size=3) + xlab("Minutes") + ylab("PI (Events)")
  }
  show(p)
  
  if(SaveToFile==TRUE){ 
    graphics.off()
  }
}


##### Data output ######
## Current Type options are (case sensitive): BaselinedData, Durations, TimeBtw, TotalFeeding
OutputData.Monitors<-function(monitors,parameters,expDesign=NA,range=c(0,0),Type="BaselinedData",filename=NA){
  if(Type=="BaselinedData"){
    if(is.na(filename)){
      OutputBaselinedData.Monitors(monitors,parameters,range)
    }
    else{
      OutputBaselinedData.Monitors(monitors,parameters,range,filename)
    }
  }
  else if(Type=="Durations"){
    if(is.na(filename)){
      OutputDurationData.Monitors(monitors,parameters,expDesign,range)
    }
    else{
      OutputDurationData.Monitors(monitors,parameters,expDesign,range,filename)
    }
  }
  else if(Type=="TimeBtw"){
    if(is.na(filename)){
      OutputIntervalData.Monitors(monitors,parameters,expDesign,range) 
    }
    else{
      OutputIntervalData.Monitors(monitors,parameters,expDesign,range,filename) 
    }
  }
  else if(Type=="TotalFeeding"){
    if(is.na(filename)){
      OutputTotalFeeding.Monitors(monitors,parameters,expDesign,range)
    }
    else {
      OutputTotalFeeding.Monitors(monitors,parameters,expDesign,range,filename)
    }
  }
  else {
    stop("Data output type does not exist.")
  }
}

##### Data functions for Exploring Individual DFMs ######
Feeding.Summary.DFM<-function(dfm,range=c(0,0),TransformLicks=TRUE){
  if(dfm$Parameters$Chamber.Size==1)
    Feeding.Summary.OneWell(dfm,range,TransformLicks)
  else if(dfm$Parameters$Chamber.Size==2)
    Feeding.Summary.TwoWell(dfm,range,TransformLicks)
  else
    stop("Feeding Summary not implemented for this DFM type.")    
}
BinnedFeeding.Summary.DFM<-function(dfm,binsize.min,range=c(0,0),TransformLicks=TRUE){
  if(sum(range)!=0) {
    m.min<-range[1]
    m.max<-range[2]
  }    
  else{
    m.min<-0
    m.max<-max(dfm$RawData$Minutes)
  }
  if(m.min>m.max)
    m.max<-m.min+1
  
  y<-seq(m.min,m.max,by=binsize.min)
  if(y[length(y)]<m.max)
    y<-c(y,m.max)
  
  tmpMatrix<-cbind(y[-length(y)],y[-1])
  intervals<-cut(y+0.000001,y,include.lowest=TRUE,dig.lab=8)
  intervals<-intervals[-length(intervals)]
  result<-Feeding.Summary.DFM(dfm,range=tmpMatrix[1,],TransformLicks)
  Interval<-rep(intervals[1],nrow(result))
  Minutes<-rep(mean(tmpMatrix[1,]),nrow(result))
  result<-data.frame(Interval,Minutes,result)
  for(i in 2:nrow(tmpMatrix)){
    tmp<-Feeding.Summary.DFM(dfm,range=tmpMatrix[i,],TransformLicks)
    Interval<-rep(intervals[i],nrow(tmp))
    Minutes<-rep(mean(tmpMatrix[i,]),nrow(tmp))
    tmp<-data.frame(Interval,Minutes,tmp)
    result<-rbind(result,tmp)
  }
  result
}
GetIntervalData.DFM<-function(dfm,range=c(0,0)){
  for(i in 1:12){
    tmp<-GetIntervalData.Well(dfm,i,range)
    if(i==1)
      result<-tmp
    else
      result<-rbind(result,tmp)  
  }
  result
}
GetDurationData.DFM<-function(dfm,range=c(0,0)){
  for(i in 1:12){
    tmp<-GetDurationData.Well(dfm,i,range)
    if(i==1)
      result<-tmp
    else
      result<-rbind(result,tmp)  
  }
  result
}

##### Plots for Exploring Individual DFMs ######
BinnedLicksPlot.DFM<-function(dfm,binsize.min=30,range=c(0,0),TransformLicks=TRUE){
  if(dfm$Parameters$Chamber.Size==1)
    PlotBins.Licks.DFM.OneWell(dfm,binsize.min,range,TransformLicks)
  else if(dfm$Parameters$Chamber.Size==2)
    PlotBins.Licks.DFM.TwoWell(dfm,binsize.min,range,TransformLicks)
  else
    stop("Binned lick plots not implemented for this DFM type.")    
}
BinnedDurationsPlot.DFM<-function(dfm,binsize.min=30,range=c(0,0)){
  if(dfm$Parameters$Chamber.Size==1)
    PlotBins.Durations.DFM.OneWell(dfm,binsize.min,range)
  else if(dfm$Parameters$Chamber.Size==2)
    PlotBins.Durations.DFM.TwoWell(dfm,binsize.min,range)
  else
    stop("Binned durations plots not implemented for this DFM type.")    
}
CumulativeLicksPlots.DFM<-function(dfm,SinglePlot=FALSE,TransformLicks=TRUE){
  tmp<-Feeding.Durations.Well(dfm,1)
  if(length(tmp)==1 && tmp[1]==0){
    fake.entry1<-c(FirstSampleData(dfm)$Minutes,0,0,0,0,0,0,0)
    fake.entry2<-c(LastSampleData(dfm)$Minutes,0,0,0,0,0,0,0)
    fake.entry<-data.frame(rbind(fake.entry1,fake.entry2),row.names=c(1,2))
    names(fake.entry)<-c("Minutes","Licks","Duration","TotalIntensity","AvgIntensity","MinIntensity","MaxIntensity","VarIntensity")
    tmp<-fake.entry
    
  }
  SumLicks<-cumsum(tmp$Licks)
  if(TransformLicks==TRUE)
    SumLicks<-SumLicks^0.25
  row<-1
  col<-1
  Row<-rep(row,nrow(tmp))
  Col<-rep(col,nrow(tmp))
  Well<-rep(1,nrow(tmp))
  tmp<-data.frame(tmp,SumLicks,Row,Col,Well)
  for(i in 2:12){
    if(i%%2==1) row<-row+1
    col<-col+1
    if(col>2) col<-1
    tmp2<-Feeding.Durations.Well(dfm,i)
    if(is.data.frame(tmp2)) {
      x<-tmp2$Minutes
      y<-cumsum(tmp2$Licks)
      SumLicks<-y
      if(TransformLicks==TRUE)
        SumLicks<-SumLicks^0.25
      Row<-rep(row,length(x))
      Col<-rep(col,length(x))
      Well<-rep(i,length(x))
      tmp2<-data.frame(tmp2,SumLicks,Row,Col,Well)
      tmp<-rbind(tmp,tmp2) 
    }
    else {
      fake.entry1<-c(FirstSampleData(dfm)$Minutes,0,0,0,0,0,0,0,0,row,col,i)
      fake.entry2<-c(LastSampleData(dfm)$Minutes,0,0,0,0,0,0,0,0,row,col,i)
      fake.entry<-data.frame(rbind(fake.entry1,fake.entry2))
      names(fake.entry)<-names(tmp)
      tmp<-rbind(tmp,fake.entry) 
    }
  }
  if(dfm$Parameters$Chamber.Size==2){
    
    even<-(as.numeric(as.character(tmp$Well))%%2==0)
    if(dfm$Parameters$PI.Multiplier==1){
      WellTC<-rep("WellA",nrow(tmp))
      WellTC[even]<-"WellB"
    }
    else {
      WellTC<-rep("WellB",nrow(tmp))
      WellTC[even]<-"WellA"
    }
    tmp<-data.frame(tmp,WellTC)
  }
  
  if(TransformLicks==TRUE)
    ylabel="Transformed Cumulative Licks"
  else
    ylabel="Cumulative Licks"
  if(SinglePlot==FALSE) {
    if(dfm$Parameters$Chamber.Size==1) {
      gp<-ggplot(tmp,aes(Minutes,SumLicks,color=factor(Well))) + geom_line() + facet_grid(rows=vars(Row),cols=vars(Col)) +geom_point() +
        ggtitle(paste("DFM",dfm$ID)) + ylab(ylabel) + labs(color="Chamber")
    }
    else {
      gp<-ggplot(tmp,aes(Minutes,SumLicks,color=factor(Well))) + geom_line() + facet_grid(rows=vars(Row),cols=vars(WellTC)) +geom_point() +
        ggtitle(paste("DFM",dfm$ID)) + ylab(ylabel) + labs(color="Chamber")
    }
    
  }
  else {
    gp<-ggplot(tmp,aes(Minutes,SumLicks,color=factor(Well))) + geom_line(size=1.2) +
      ggtitle(paste("DFM",dfm$ID))+ ylab(ylabel)+ labs(color="Chamber")
  }
  show(gp)
}


CumulativePI.DFM<-function(dfm, range=c(0,0), ShowPlots=TRUE, SinglePlot=FALSE){
  ## Get the Feeding.PI
  chambers<-1:nrow(dfm$Parameters$Chamber.Sets)
  Minutes<-Minutes(dfm,range)
  for(i in 1:nrow(dfm$Parameters$Chamber.Sets)) {
    
    if(dfm$Parameters$PI.Multiplier==1){
      wellA<-dfm$Parameters$Chamber.Sets[i,1]
      wellB<-dfm$Parameters$Chamber.Sets[i,2] 
    }
    else {
      wellB<-dfm$Parameters$Chamber.Sets[i,1]
      wellA<-dfm$Parameters$Chamber.Sets[i,2] 
    }
    
    a<-FeedingData.Well.Licks(dfm,wellA,range)
    b<-FeedingData.Well.Licks(dfm,wellB,range)
    FeedingLicksA<-cumsum(a)
    FeedingLicksB<-cumsum(b)
    
    ## Here it is the instantaneous PI
    Feeding.PI<-data.frame(Minutes,(FeedingLicksA - FeedingLicksB)/(FeedingLicksA+FeedingLicksB),FeedingLicksA+FeedingLicksB)
    Feeding.PI<-Feeding.PI[a+b>0,]
    Feeding.PI<-data.frame(Feeding.PI,rep(chambers[i],nrow(Feeding.PI)))
    names(Feeding.PI)<-c("Minutes","PI","Licks","Chamber")
    
    if(i==1)
      results<-Feeding.PI
    else
      results<-rbind(results,Feeding.PI)
    
  }
  if(ShowPlots==TRUE){
    if(SinglePlot==FALSE) 
      gp<-ggplot(results,aes(Minutes,PI,color=Licks)) + geom_line() + facet_grid(rows=vars(factor(Chamber))) +geom_point() +
        ggtitle(paste("DFM",dfm$ID)) + ylab("PI (Licks)") + labs(color="Licks") + ylim(c(-1,1)) + scale_color_gradientn(colours = rainbow(5))
    else
      gp<-ggplot(results,aes(Minutes,PI,color=factor(Chamber))) + geom_line() +geom_point() +
        ggtitle(paste("DFM",dfm$ID)) + ylab("PI (Licks)") + labs(color="Chamber") + ylim(c(-1,1))
    
    show(gp)
  }
  results
}
CumulativeEventPI.DFM<-function(dfm, events.limit=NA, range=c(0,0), ShowPlots=TRUE, SinglePlot=FALSE){
  ## Get the Feeding.PI
  chambers<-1:nrow(dfm$Parameters$Chamber.Sets)
  Minutes<-Minutes(dfm,range)
  for(i in 1:nrow(dfm$Parameters$Chamber.Sets)) {
    
    if(dfm$Parameters$PI.Multiplier==1){
      wellA<-dfm$Parameters$Chamber.Sets[i,1]
      wellB<-dfm$Parameters$Chamber.Sets[i,2] 
    }
    else {
      wellB<-dfm$Parameters$Chamber.Sets[i,1]
      wellA<-dfm$Parameters$Chamber.Sets[i,2] 
    }
    
    a<-FeedingData.Well.Events(dfm,wellA,range)
    b<-FeedingData.Well.Events(dfm,wellB,range)
    a[a>0]<-1
    b[b>0]<-1
    FeedingEventsA<-cumsum(a)
    FeedingEventsB<-cumsum(b)
    
    ## Here it is the instantaneous PI
    Feeding.PI<-data.frame(Minutes,(FeedingEventsA - FeedingEventsB)/(FeedingEventsA+FeedingEventsB))
    Feeding.PI<-Feeding.PI[a+b>0,]
    if(nrow(Feeding.PI)>0){
      Feeding.PI<-data.frame(Feeding.PI,rep(chambers[i],nrow(Feeding.PI)),1:nrow(Feeding.PI))
      names(Feeding.PI)<-c("Minutes","PI","Chamber","EventNum")
      if(is.na(events.limit)==FALSE && events.limit<nrow(Feeding.PI)){
        Feeding.PI<-Feeding.PI[1:events.limit,]
      }
      if(i==1)
        results<-Feeding.PI
      else
        results<-rbind(results,Feeding.PI)
    }
  }
  if(ShowPlots==TRUE){
    if(SinglePlot==FALSE) 
      gp<-ggplot(results,aes(Minutes,PI,color=EventNum)) + geom_line() + facet_grid(rows=vars(factor(Chamber))) +geom_point() +
        ggtitle(paste("DFM",dfm$ID)) + ylab("PI (Events)") + labs(color="Events") + ylim(c(-1,1))+ scale_color_gradientn(colours = rainbow(5))
    else
      gp<-ggplot(results,aes(Minutes,PI,color=factor(Chamber))) + geom_line() + geom_point() +
        ggtitle(paste("DFM",dfm$ID)) + ylab("PI (Events)") + labs(color="Events") + ylim(c(-1,1))
    show(gp)
  }
  results
}



