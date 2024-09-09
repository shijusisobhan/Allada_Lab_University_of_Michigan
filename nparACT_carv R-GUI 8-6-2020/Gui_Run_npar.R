



#****************************************************************************************************

# Tis is the function for run nparact and carv analysis for subject and return the results to server
# usage eg:-

#st.day<- '25/12/2015'
#ed.day<- '29/12/2015'
#st.time <- '12:00:00 PM'
# ed.time<- '12:00:00 PM'
# Ts<-1
#PATH <- '/Desktop/myfolder'
#Subject.ID <- '1741-001'


# Gui_Run_npar(st.day, ed.day, st.time, ed.time, Ts, PATH, Subject.ID)

#*****************************************************************************************************

Gui_Run_npar<-function(st.day, ed.day, st.time, ed.time, Ts, PATH, Subject.ID){
  
  
  
  st.day<<-st.day
  
  ed.day<<-ed.day
  
  st.time<<-st.time
  
  ed.time<<-ed.time
  
  Ts<<-Ts
  
  subject.ID<<-Subject.ID
  
  PATH<<-PATH
  
  
  #*********************************************************************************
  
  Result.info<-data.frame(Subject.ID, st.day,st.time, ed.day, ed.time, Ts)
  
  
  colnames(Result.info)<-c("subject", "start.day","start.time", "end.day",
                              "end.time", "Sampling.time ")
  
  
#**********************************************************************************
  

source("nparACT_carv.R")


# sub data for particular subject
actiware.data.Merge_sub<- subset(actiware.data.Merge, actiware.data.Merge$sub.ID==Subject.ID)
actiware.data<<- readActiwareCSV(actiware.data.Merge_sub)



# *************************nparACT analysis **************************************************

SR=1/(as.numeric(Ts)*60) # satmpling rate

if (SR<(1/60)){
  
sleepstudy.activity<<-actiware.data[,c("DayTime", "Activity")]
sleepstudy.WL<<-actiware.data[,c("DayTime", "White.Light")]


sleepstudy.activity<<-na.omit(sleepstudy.activity)
sleepstudy.WL<<-na.omit(sleepstudy.WL)


source("npar_base_New.R") # This is the modified version of npar_base function

npar_plot_data<<-npar_base_New("sleepstudy.activity", SR, cutoff = 1)
Result_npar_Activity<<-npar_plot_data[[4]]

Result_npar_Activity$IS<-NaN
Result_npar_Activity$IV<-NaN
Result_npar_Activity$RA<-NaN
Result_npar_Activity$L5<-NaN
Result_npar_Activity$L5_starttime<-NaN
Result_npar_Activity$M10<-NaN
Result_npar_Activity$M10_starttime<-NaN
npar_plot_data[[1]]=NaN


npar_plot_data_WL<<-npar_base_New("sleepstudy.WL", SR, cutoff = 1)
Result_npar_WL<<-npar_plot_data_WL[[4]]


Result_npar_WL$IS<-NaN
Result_npar_WL$IV<-NaN
Result_npar_WL$RA<-NaN
Result_npar_WL$L5<-NaN
Result_npar_WL$L5_starttime<-NaN
Result_npar_WL$M10<-NaN
Result_npar_WL$M10_starttime<-NaN


}



else{
sleepstudy.activity<<-actiware.data[,c("DayTime", "Activity")]
sleepstudy.WL<<-actiware.data[,c("DayTime", "White.Light")]




sleepstudy.activity<<-na.omit(sleepstudy.activity)
sleepstudy.WL<<-na.omit(sleepstudy.WL)


source("npar_base_New.R") # This is the modified version of npar_base function

npar_plot_data<<-npar_base_New("sleepstudy.activity", SR, cutoff = 1)
Result_npar_Activity<<-npar_plot_data[[4]]





npar_plot_data_WL<<-npar_base_New("sleepstudy.WL", SR, cutoff = 1)
Result_npar_WL<<-npar_plot_data_WL[[4]]

}



colnames(Result_npar_Activity)<-c("Activity.IS","Activity.IV", "Activity.RA", 
                                  "Activity.L5", "Activity.L5ST", 
                                  "Activity.M10", "Activity.M10ST")

colnames(Result_npar_WL)<-c("WLight.IS","WLight.IV", "WLight.RA", 
                                  "WLight.L5", "WLight.L5ST", 
                                  "WLight.M10", "WLight.M10ST")


# *************************carv analysis **************************************************



my.linearFit <<- lActReg(actiware.data) ## lACTREG ()  called in fitNLS

my.nonlinearFit <<- fitNLS(actiware.data)

carv.Result<<-NLSredux(actiware.data)

# NLSredux(actiware.data)

ID="AB"



# Retuning the results including npar, carv and data for ploting

return(list(Result.info,Result_npar_Activity,Result_npar_WL,carv.Result,ID,actiware.data,my.nonlinearFit,
            
            npar_plot_data[[1]], npar_plot_data[[2]],npar_plot_data[[3]], SR))

}


