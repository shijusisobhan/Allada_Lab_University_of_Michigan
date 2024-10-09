# Load necessary library


rm(list=ls())
library(dplyr)


 setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/SleepMat_all/YK/FLIC_data_Eduction/FLIC_2_Sleep/YK_Flic2Sleep_Data')

# Read the .txt file into R (assuming tab-separated, adjust as needed)
data <- read.table("Monitor3.txt", header = F, sep = "\t")



start_date <- as.character(head(data[2],1))
end_date <- as.character(tail(data[2],1))
start_time <- as.character(head(data[3],1))
end_time <- as.character(tail(data[3],1))

# Combine the start date and time, and end date and time into POSIXct objects
start_datetime <- as.POSIXct(paste(start_date, start_time), format="%d %b %Y %H:%M:%S")
end_datetime <- as.POSIXct(paste(end_date, end_time), format="%d %b %Y %H:%M:%S")

second(start_datetime) <- 0
second(end_datetime) <- 0

# Generate a sequence of times at one-minute intervals
time_sequence <- seq(from = start_datetime, to = end_datetime, by = "1 min")

# Create a data frame with two columns: Date and Time
result_df <- data.frame(
  Date = format(time_sequence, "%d %b %Y"),  # Format date as "13 May 2024"
  Time = format(time_sequence, "%H:%M:%S")            # Extract the time
)

data[3]<- result_df$Time

write.table(data, 'Monitor3.txt', sep = "\t", row.names = F, col.names = F, quote = FALSE)