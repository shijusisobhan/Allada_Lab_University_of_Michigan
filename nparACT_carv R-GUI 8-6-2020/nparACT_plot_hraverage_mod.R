
#********************************************************************************
# Program for Actigraphy plot hour-wise average
# Eg:-nparACT_plot_hraverage_mod(data, minaverage, start.time, a, SR)
#*****************************************************************************


nparACT_plot_hraverage_mod = function(data, minaverage, start.time, a, SR){
  hraverage <- matrix(NA)
  for (i in 1:24){
    hraverage [i] <- mean(minaverage[(((i-1)*60)+1):(60*i)])
  }
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
  start.time <- as.numeric(timeinfo[1])  
  
  for_hraverage_plot.time <- rep(seq(1,24),2)
  seq <- seq(start.time, length.out = 24)
  if (seq[1] == 0){
    seq[1] = 24 ## changed 21 12 16, accounts for case in which start time is midnight
  }
  hraverage_plot_time <- for_hraverage_plot.time[seq]
  hraverage_plot_time[hraverage_plot_time==24] <- 0
  hraverage_plot_df <- data.frame (hraverage_plot_time, hraverage)
  
  ppp <- ggplot(hraverage_plot_df, aes(x=hraverage_plot_time, y = hraverage))+ 
    geom_bar(stat="identity", width = 1, position = position_dodge(width = 0.5))+
    theme_bw()+
    scale_x_discrete(expand = c(0,0), limits = c(seq(0:23)-1), breaks = seq(0,23))+
    expand_limits(x=-1)+
    xlab("Time \n (Start: 0am)")+
    ylab("Movement Intensity")+
    ggtitle(paste("Hour-wise average across days"))+
    #ggtitle("Actigraphy Plot (24 hrs)", subtitle = "Hour-wise average across days")+
    theme(plot.margin=unit(c(1,2,1,2),"cm"),
          axis.text.x = element_blank(),
          axis.title.x = element_text(face="bold", size=12),
          axis.title.y = element_text(face="bold", size=12, vjust = 1.2),
          plot.title = element_text(face="bold", size=14, vjust = 1))
  
}