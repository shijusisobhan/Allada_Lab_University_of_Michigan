

setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/SleepMat_all/YK/FLIC_data_Eduction/FLIC_2_Sleep/YK_Flic2Sleep_Data/DAM1234')

library(readxl)
library(dplyr)
library(openxlsx)
library(gridExtra)


sleepdata<- read_excel('DAM1234_Sleep_interval.xls')

Activity_data<- read_excel('DAM1234_Activity_interval.xls')


Activated<-sleepdata[25,-c(1:5)]
control<-sleepdata[55,-c(1:5)]
time_data<-seq(0.5,48,0.5)



plot(time_data, Activated, type="o",col="red", xlab='ZT', ylab='sleep/30 min')
lines(time_data,control, 'black')
segments(x0 = 24, y0 = 1, x1 = 36, y1 = 1, col = "blue", lwd = 5)




Activated<-Activity_data[25,-c(1:5)]
control<-Activity_data[55,-c(1:5)]
time_data<-seq(0.5,48,0.5)

p2<-plot(time_data, Activated, type="o",col="red", xlab='ZT', ylab='Activity/30 min')
lines(time_data,control, 'black')
segments(x0 = 24, y0 = 1, x1 = 36, y1 = 1, col = "blue", lwd = 5)

