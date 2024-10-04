rm(list=ls())
#setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/SleepMat_all/YK/FLIC_data_Eduction')
setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/SleepMat_all/YK/FLIC_data_Eduction/FLIC_test')

library(ggplot2)
library(stats)
library(gridExtra)
library(reshape2)
library(gtools)

source("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/GitHub_folder/Allada_Lab_University_of_Michigan/FlIC_2_SLEEP/FLIC_data_analysis_SW.R")

DFM_number<-4

p1<-ParametersClass.SingleWell(Feeding.Threshold = 10)
expDesign<-read.csv("ExpDesign.csv")
Binned_data<-BinnedFeeding.Summary.Monitors(DFM_number,p1,1,expDesign = expDesign, TransformLicks=FALSE)
flic_data<-Binned_data$Results


library(tidyverse)

flic_data<-flic_data[,c(3,5:6)]

flic_data_wide<-pivot_wider(flic_data,
                            names_prefix = 'Ch',
                            names_from = 'Chamber',
                            values_from = 'Licks')





start_date <- head(DFM4$RawData$Date,1)
end_date <- tail(DFM4$RawData$Date,1)
start_time <- head(DFM4$RawData$Time,1)
end_time <- tail(DFM4$RawData$Time,1)

# Combine the start date and time, and end date and time into POSIXct objects
start_datetime <- as.POSIXct(paste(start_date, start_time), format="%m/%d/%Y %H:%M:%S")
end_datetime <- as.POSIXct(paste(end_date, end_time), format="%m/%d/%Y %H:%M:%S")

# Generate a sequence of times at one-minute intervals
time_sequence <- seq(from = start_datetime, to = end_datetime, by = "1 min")

# Create a data frame with two columns: Date and Time
result_df <- data.frame(
  Date = format(time_sequence, "%d %b %Y"),  # Format date as "13 May 2024"
  Time = format(time_sequence, "%H:%M:%S")            # Extract the time
)


Monitor_data<-cbind(result_df, flic_data_wide)
Monitor_data[paste0('Ch',13:32)]<-0


# Save the data frame as a .txt file
write.table(Monitor_data, file = paste0('Monitor',DFM_number,'.txt'), sep = "\t", row.names = T, col.names = F, quote = FALSE)



