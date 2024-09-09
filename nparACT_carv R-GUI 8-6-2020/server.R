
#*******************************************************************************************

# (c) Shiju Sisobhan (shiju.sisobhan@northwestern.edu)

# shiny server strts
# It have different element (functions)

   # 1. reactive- It is a fuction that perform an operation
        # It can be called in the server by using its name Eg; Reactive1()

   # 2. outtput--> It return the outputs to the GUI. It have different structure
         # a)renderText<- return text output
         # b)renderTable<--return table ouput
         # c)renderPlot<- return plot output

# 3. downloadHandler<-- It is used for downloading data and plots in either csv or pdf format
     

#**********************************************************************************************


shinyServer(function(input, output){
  
  

  # Find out all the subject in the uploaded file
  
  output$subject.N<- renderTable({
    
   files_load<-input$file$datapath # load the actigraphy files


    if(is.null(files_load)){return(NULL)}
     else{NM1= subject_finder(files_load)} # This return the subject present in the uploaded files
   
   
    
  })
  
#*********************************************************************
  

  
  
  # specification file output and disable manual inputs
  
  
  output$Run.spec<-renderTable({
    
    
    sp.file <- input$spec.file$datapath # load specification file (csv)

    if(is.null(sp.file)){
      {enable("subject")}
      return(NULL)}
    else
    {disable("subject")
      disable("start.day")
      disable("end.day")
      disable("start.time")
      disable("end.time")}
    Run_spec<-read.csv(sp.file) # Give specification for analysis (subject, strt date & time.......)
   

    
    
    
  })
  
  
  #*********************************************************************
  
   #Reactive expression for shiny start
   # Reactive element can be called inside the server by their name
   # eg: Subject.files.merging()
  
  #********************************************************************
  

  
  # Merge all the actigraphy files 
  
  Subject.files.merging<-reactive({
    
    
    PATH<- input$file$datapath # Actgraphy files
    
    if(is.null(PATH)){return(NULL)}
    
    else{
    my.CSV.files<-PATH

    actiware.data.Merge <<- readActiwareCSVs_action(my.CSV.files) # It returns merged version of all actigraphy files

    actiware.data.Merge$time <- .timeString24(actiware.data.Merge$timeString) # Add an additional colum

    actiware.data.Merge}
    
    
  })
  
  # eg: tmp<-Subject.files.merging()
  
  
  
  
  # This returns all the information of subject (start& end time, sampling time.......)
  
  Run_subject_info <- reactive({
  
 
    
   sp.file <- input$spec.file$datapath
    
    
    if (is.null(sp.file)) {
      
      Data.merge=Subject.files.merging()
      subject.ID<- input$subject
      
      Run_sub=list()
      Run_sub[[1]]=Sujcet_info_finder(Data.merge,subject.ID)
      Result_sub<-Run_sub
      
      
    } else {
      
      Data.merge=Subject.files.merging()
      
      spec.info<-read.csv(sp.file)
      
      N.sub<<-nrow(spec.info) # Multisubject
      
      
      Run_sub=list()
      
      
      for (i in 1:N.sub){
        
        
        subject.ID<- as.character(spec.info$subject[i])
        
        Run_sub[[i]]=Sujcet_info_finder(Data.merge, subject.ID)
        
      }
      
      Result_sub<<- do.call("rbind",Run_sub)
      
      
      
    }
 
    
  })
  
 
  #  Eg: Run_subject_info() 
  
  
  
 # This returns an error message if the subject is not present in the uploded files
  
 error.msg<-reactive({
   
   R1<-Run_subject_info()
   
   if('Error' %in% Result_sub){return('Error! Subject not Uploaded')
     
   }
   
 }) 
  
  
  
#  Error message output in the ui 
 
output$txt_error<-renderText({
  
 input$action

  
  tryCatch(isolate(paste("", error.msg())),
           
           error=function(e){return(NULL)})
   
  
})



  
# subject information output

output$subject.information<-renderTable({
      
      input$action

      tryCatch(isolate(Run_subject_info()),

                error=function(e){return(NULL)})
  
      
      
    })
  
  
  
  
  


#**********************************************************************************
  
  # Do the carv & nparACT analysis for the subject

# This return carv, nparact output
  
  Run_npar1 <- reactive({
    
    
    sp.file <- input$spec.file$datapath

    
    
    if (is.null(sp.file)) {
      
      subject.ID<<- input$subject
      Run_npar=list()
      N.sub<<-1 # One subject
      st.day <- input$start.day
      ed.day <- input$end.day
      st.time <- input$start.time
      ed.time <- input$end.time
      
      
      Ts <- input$sampling.time
      
      
      PATH<- input$file$datapath
      
      Run_npar[[1]]=Gui_Run_npar(st.day, ed.day, st.time, ed.time, Ts, PATH,  subject.ID) 
      
      Result_all<<-Run_npar
      
    } else {
      
      spec.info<-read.csv(sp.file)
      
      N.sub<<-nrow(spec.info) # Multisubject
     
      
      Run_npar=list()
      
      # loop for multi subject analysis
      
      for (i in 1:N.sub){
      
      st.day <- as.character(spec.info$start.day[i])
      ed.day <- as.character(spec.info$end.day[i])
      st.time <- as.character(spec.info$start.time[i])
      ed.time <- as.character(spec.info$end.time[i])
      subject.ID<<- as.character(spec.info$subject[i])
      Ts <- input$sampling.time
      
      
      PATH<- input$file$datapath
      
      Run_npar[[i]]=Gui_Run_npar(st.day, ed.day, st.time, ed.time, Ts, PATH,  subject.ID)
      
      }
      
     Result_all<<-Run_npar
      
    }
    
    

    
  })
  
  
  
  # Eg: Run_npar1()
  
  
  
  
  ## ****************************************************************
  # carv & nparact analysis for Individual file
  
  Run_Individual <- reactive({
    
    PATH <- input$file$datapath
    File_names<-input$file$name
    
    if(is.null(PATH)){return(NULL)}
    else{
    Run_Individual<-Run_carv_individual(PATH, File_names, input$ind_st_T)}
    
  }) 
  
  
  
Run_carv <- reactive({
  

  
  Run_carv<-Run_Individual()
  
  Run_carv[1]
  
})


Run_npar <- reactive({
  
  PATH <- input$file$datapath
  File_names<-input$file$name
  
  Run_npar<-Run_Individual()
  
  Run_npar[2]
  
})


# Display the individual file results

output$Individual_Result_carv<-renderTable({
  
  input$I_Run
  
  isolate(Run_carv())
  
})

# Eg: Run_Individual()




 # Individual file nparACT results

output$Individual_Result_napar<-renderTable({
  
  input$I_Run
  
  isolate(Run_npar())
  
})



##*********************************************************************************

 # Display nparACT and carv Output for each subject
  
# nparact for activity
  
  output$Run_npar_Activity <- renderTable({
    
    tryCatch({NM1= Run_npar1()
    Run_npar_Activity=list()
    for (i in 1:N.sub) {

      Run_npar_Activity[[i]]=NM1[[i]][[2]]
    }

    df1 <- do.call("rbind",Run_npar_Activity)},

    error=function(e){return()})
    
    
   
  })
  
  
  
 # warning message if the TS>1 
  
  warning.msg<-reactive({
    

    
    if(as.numeric(input$sampling.time)>1){return('Warning! nparACT not works if sampling time>1')
      
    }
    
  })
  
  
  # Display warning message
  
  output$txt_warning<-renderText({
  

    tryCatch(paste("", warning.msg()),
             
             error=function(e){return(NULL)})
    
    
  }) 
  
  
  
  #  nparACT for WL
  
  output$Run_npar_WL <- renderTable({
    
    tryCatch({NM1= Run_npar1()
    Run_npar_WL=list()
    for (i in 1:N.sub) {

      Run_npar_WL[[i]]=NM1[[i]][[3]]
    }

    df1 <- do.call("rbind",Run_npar_WL)},



    error=function(e){return(NULL)}
    )
    


    
    
    
  })
  
  
  
  # carv output
  
  output$Result_carv<- renderTable({
    

    tryCatch({NM1= Run_npar1()
    Run_carv=list()
    for (i in 1:N.sub) {

      Run_carv[[i]]=NM1[[i]][[4]]

    }

    df1 <- do.call("rbind",Run_carv)},

    error=function(e){return(NULL)}

    )

    
 
    
  })
  
  
  
  ##*****************************************************************************
  
  # Diwnload the output files

  
  output$Result.csv <- downloadHandler(
      filename = function() {
        "Result.csv"
      },
      content = function(file) {
        Run_info=list()
        Run_npar_Activity=list()
        Run_npar_WL=list()
        Run_npar_carv=list()
        
        for (i in 1:N.sub) {
          Run_info[[i]]=Result_all[[i]][[1]]
          Run_npar_Activity[[i]]=Result_all[[i]][[2]]
          Run_npar_WL[[i]]=Result_all[[i]][[3]]
          Run_npar_carv[[i]]=Result_all[[i]][[4]]
          
        }
        
        
        df1 <- do.call("rbind",Run_info)
        df2 <- do.call("rbind",Run_npar_Activity)
        df3 <- do.call("rbind",Run_npar_WL)
        df4 <- do.call("rbind",Run_npar_carv)
        
        out_all<-cbind(df1,df4,df2,df3)
        
        write.csv(out_all, file, row.names = F, sep=',')
        
       
        

      }
    )
  
  
  
  
  output$Result_individual.csv <- downloadHandler(
    filename = function() {
      "Result_individual.csv"
    },
    content = function(file) {
      NM2= Run_Individual()
      
     df<-NM2
      
      write.csv(df, file, row.names = F, sep=',')
      
      
      
    }
  )
  
  
  
## *************************************************************************************************

    
# Download the plots
  

          output$plot.pdf = downloadHandler(
            filename = 'plot.pdf',
            content = function(file) {
              pdf(file)

              
              
              
              for (i_Nsub in 1:N.sub){
                

                id="AB"
                actidat11<<-Result_all[[i_Nsub ]][[6]]
                nlsfit11<<-Result_all[[i_Nsub ]][[7]]
                
                
                
                plotDF11<- actidat11[actidat11$idcode==id,c("time","cloktime","lActivity", "Time", "dateString")]
                plotDF11 <- na.omit(plotDF11)
                
                plotDF11$plottime<-(plotDF11$cloktime-plotDF11$time[1])%%24
                
     
                plotDF11$cloktime<-(plotDF11$cloktime-plotDF11$time[1])
                
                plotDF11$day <- (plotDF11$cloktime%/%24)+1
                plotDF11$predlact <- predict(nlsfit11[[id]],plotDF11$cloktime)
                
                
                plotDF111<<-plotDF11
                
                plotDF1 <<- split(plotDF11,plotDF11$day) # split data for each day
                

                
                plist <- list()
                

                id="AB"
                
                max_plot<-length(plotDF1)
                
                type="p"
                t1=0
                t2=24
                
                plotlabel=c()
                
                
                test.data1<<-plotDF1[[1]]
                
                
                for (ilab in 1:5){
                  plotlabel[ilab] <- paste(round((plotDF1[[1]]$time[1]+(ilab-1)*6)%%24))
                }
                
                
                
                for (i in 1:max_plot){
                  
                  
                  day <- i
                  

 
                  plist[[day]] <- ggplot(data = plotDF1[[day]]) + 
                    geom_point(aes(x = plottime, y = lActivity)) +
                    geom_line(aes(x = plottime, y = predlact), col=2,lwd=2) +
                    scale_x_continuous(name = "Clock time", limits = c(0,24), breaks=c(0,6,12,18,24), labels =plotlabel)+
                    ggtitle(paste("Sub.ID", Result_all[[i_Nsub]][[1]]$subject,"," , "  day-",plotDF1[[day]]$dateString[1]))+
                    xlab("Clock time")+
                    theme_bw()
                  
 
                  
                  
                }
                
 
                
                pl<- lapply(1:max_plot, function(i){
                  
                  (plist[[i]])
                  
                })
                
                
                
                sp.file <- input$spec.file$datapath
                
                
                
                if (is.null(sp.file)) {
                  
                  ml <- marrangeGrob(pl, nrow=4, ncol=2,
                                     top = paste("Date & Time:",st.day,
                                                 " ",st.time, "-", ed.day, " ",
                                                 ed.time, " Ts=", input$sampling.time, " min"))
                  
                } else {
                  
                  spec.info<-read.csv(sp.file)
                  
                                
                    
                    st.day <- as.character(spec.info$start.day[i_Nsub])
                    ed.day <- as.character(spec.info$end.day[i_Nsub])
                    st.time <- as.character(spec.info$start.time[i_Nsub])
                    ed.time <- as.character(spec.info$end.time[i_Nsub])
                    
                    
                    
                    ml <- marrangeGrob(pl, nrow=4, ncol=2,
                                       top = paste("Date & Time:",st.day,
                                                   " ",st.time, "-", ed.day, " ",
                                                   ed.time, " Ts=", input$sampling.time, " min"))
                    
                  }
                
                
          

            print(ml)
                
                
                
       # nparACT plot
            
           start.time.plot=NULL
           
           tryCatch({
             
             ml1=list()
             
             
             ml1[[1]]<-nparACT_plot_minaverage_mod(Result_all[[i_Nsub ]][[8]], 
                                                   Result_all[[i_Nsub ]][[9]],
                                                   start.time.plot ,Result_all[[i_Nsub ]][[10]],
                                                   Result_all[[i_Nsub ]][[11]], 
                                                   Result_all[[i_Nsub]][[1]]$subject)
             
             ml1[[2]] <-nparACT_plot_hraverage_mod(Result_all[[i_Nsub ]][[8]], Result_all[[i_Nsub ]][[9]],
                                                   start.time.plot ,Result_all[[i_Nsub ]][[10]],
                                                   Result_all[[i_Nsub ]][[11]])
             
             ml2 <- marrangeGrob(ml1, nrow=2, ncol=1,
                                 top = paste("Date & Time:",st.day,
                                             " ",st.time, "-", ed.day, " ",
                                             ed.time, " Ts=", input$sampling.time, " min"))
             
             
             print(ml2) # print npar plots
             
             
           },

           error=function(e){return(NULL)}
           )
  

                
              }
            
            
            
            dev.off()

            })

          

        
})


##************************END***************************************##########