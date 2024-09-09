
#********************************************************************************
# Program for Actigraphy plot minute-wise average
# Eg:-nparACT_plot_minaverage_mod(data, minaverage, start.time, a, SR)
#*****************************************************************************


nparACT_plot_minaverage_mod = function(data, minaverage, start.time, a, SR, sub.IDI){
  daytime <- matrix(NA)
  time <- data$time
  time <- as.character(time)
  for (v in seq(1,a,(SR*60*60))){  
    daytime[v] <- time[v]
  }
  daytime <- na.omit(daytime)
  daytime <- as.character(daytime)
  
  temp = unlist(str_split(daytime, ' ') )
  temp_nums = 1:length(temp)
  timeinfo = temp[ (temp_nums %% 2) == 0 ] 
  temp = unlist(str_split(timeinfo, ':') )
  temp_nums = 1:length(temp)
  timeinfo = temp[ (temp_nums %% 3) == 1 ] 
  start.time_h <- as.numeric(timeinfo[1])
  start.time_min <- as.numeric(temp[2]) ## changed 21 12 16 -> select start minute, was [1] for hour previously
  
  for_minaverage_plot.time <- rep(seq(1,1440),2)
  seq <- seq(start.time_h*60+start.time_min, length.out = 1440) ## 21 12 16: was start.time*60 when input was hours only
  if (seq[1] == 0){  ## changed 21 12 16, accounts for case in which start time is midnight
    seq[1] = 1440
  }
  minaverage_plot_time <- for_minaverage_plot.time[seq] # time info, vectors will be joined in next step
  
  minaverage_plot_df <- data.frame (minaverage_plot_time, minaverage)
  
  pp <- ggplot(minaverage_plot_df, aes(x=minaverage_plot_time, y = minaverage))+ 
    geom_bar(stat="identity", width = 1, position = position_dodge(width = 0.5))+
    theme_bw()+
    scale_x_discrete(limits = c(seq(1:1440)), breaks = seq(1,1440,60))+
    expand_limits(x=c(-30,1470))+
    xlab("Time \n (Start: 0am)")+
    ylab("Movement Intensity")+
    ggtitle(paste("Sub.ID:", sub.IDI, ", Minute-wise average across days"))+
    #ggtitle("Actigraphy Plot (24 hrs)", subtitle = "Minute-wise average across days")+
    theme(plot.margin=unit(c(1,2,1,2),"cm"),
          axis.text.x = element_blank(),
          axis.title.x = element_text(face="bold", size=12),
          axis.title.y = element_text(face="bold", size=12, vjust = 1.2),
          plot.title = element_text(face="bold", size=14, vjust = 1))

}
