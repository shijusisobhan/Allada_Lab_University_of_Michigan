
#*****************************************************************

# This is the function returns all the information reagrding the subject

# Eg:- Sujcet_info_finder(Merged data, 1741-001')

#***************************************************************
Sujcet_info_finder<-function(data.Merge,S.ID){
  
   source("nparACT_carv.R")
    
    sub.out<-c()
    
    data.Merge.sub<- subset(data.Merge, data.Merge$sub.ID==S.ID)
    
    
 
    if (nrow(data.Merge.sub)==0){


      sub.out$Suject.ID=S.ID
      sub.out$default.Ts <- paste('Error')
      sub.out$start.date='-'
      sub.out$start.time='-'
      sub.out$end.date='-'
      sub.out$end.time='-'


      #error=function(e){stop("Subject not uploded")}


    }

else{

    sub.out$Suject.ID=S.ID
    default.sampling.time1<<-round((as.numeric(data.Merge.sub$time[2])-as.numeric(data.Merge.sub$time[1]))*60, 1)
    sub.out$default.Ts <- default.sampling.time1
    sub.out$start.date=head(data.Merge.sub$Date, n=1)
    sub.out$start.time=head(data.Merge.sub$timeString, n=1)
    sub.out$end.date=tail(data.Merge.sub$Date, n=1)
    sub.out$end.time=tail(data.Merge.sub$timeString, n=1)

}
    
    
    return(sub.out)
  
}

