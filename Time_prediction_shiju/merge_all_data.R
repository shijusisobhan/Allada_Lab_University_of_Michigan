
rm(list=ls())

setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/TimeSignature/New Project/Amplitude estimator/TimeMachine-main/Time_prediction_shiju')
library(dplyr)

DLMO_rounding<- function(Data1) {
  DLMO_round<-ceiling(floor(Data1[1,-1]*2)/2)
  DLMO_round$ext_gene<-'DLMO_Intiger'
  DLMO_round<-DLMO_round %>% select(ext_gene, everything())
  position<-1
  
  Data_new<-rbind(Data1[1:position,], DLMO_round, Data1[(position+1):nrow(Data1),])
  
}

Archer_Data<-read.csv('Archer_data_and_DLMO_all.csv')
Archer_Data_DL_round<-DLMO_rounding(Archer_Data)

Braun_data<-read.csv('Braun_data_and_DLMO_all.csv')
Braun_data_DL_round<-DLMO_rounding(Braun_data)

Mollar_Data<-read.csv('Moller_data_and_DLMO_all.csv')
Mollar_data_DL_round<-DLMO_rounding(Mollar_Data)


MAB<-Mollar_data_DL_round %>% left_join(Archer_Data_DL_round, by='ext_gene') %>% 
  left_join(Braun_data_DL_round, by='ext_gene')

write.csv(MAB, 'Moller_archer_Broun_data.csv')

