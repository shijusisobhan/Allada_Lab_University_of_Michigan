rm(list=ls())

library(dplyr)

setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/TimeSignature/New Project/Amplitude estimator/TimeMachine-main/Time_prediction_shiju')


Archer_Data<-read.csv('Archer_data_and_DLMO_all.csv')

Braun_data<-read.csv('Braun_data_and_DLMO_all.csv')

Mollar_Data<-read.csv('Moller_data_and_DLMO_all.csv')


MAB_data<-Mollar_Data %>% left_join(Archer_Data, by='ext_gene') %>% 
  left_join(Braun_data, by='ext_gene')

DLMO25<-t(MAB_data[1,-1]) # extract the DLMO value

#MAB_data<-read.csv('Moller_archer_Broun_data.csv')
MAB_data<-MAB_data[-c(1),]

JTK_results<-read.csv('JTK_z_score_all.csv')
cycling_genes<-JTK_results[which(JTK_results$BH.Q<0.05),]

TS_genes<-read.csv('TS_genes_all.csv')

TS_cyclic_genes<-cycling_genes[cycling_genes$ext_gene %in% TS_genes$ext_genes,]

MAB_TS<-MAB_data[MAB_data$ext_gene %in% TS_cyclic_genes$ext_gene,]



MAB_TS_JTK<-merge(MAB_TS, TS_cyclic_genes[c('ext_gene', 'LAG')], by.x=1 , by.y="ext_gene")

##*****************************************************************************

cbn_fun<-function(a){combn(a,2)} # Pairwise gene combination function


ratio_genes_fun<-function(a, n){
  
  B=c()
  
  for (i in 0:n-1) {B[i+1]<-a[2*i+1]-a[2*i+2]}  # (log (a/b) = log(a) - log (b))
  
  B
  
} # ratio of pairs



LAG_diff_fun<-function(a,n){
  
  B=c()
  
  for (i in 0:n-1) {
    
    B[i+1]=abs(a[2*i+1]-a[2*i+2])
    
  }
  
  B
  
} # ratio of pairs



ratio_genes_name_fun<-function(a,n){
  
  G_N=c()
  
  for (i in 0:n-1) {G_N[i+1]<-paste(a[2*i+1],',', a[2*i+2])}
  
  G_N
  
} # ratio of pairs



genes_combination_numeric<-apply(MAB_TS_JTK[,2:768],2, cbn_fun) # Pairwise combination

genes_name_combination<-cbn_fun(MAB_TS_JTK[,1])

NN<-(nrow(genes_combination_numeric))/2

#********************************************  

genes_ratio_numeric<-apply(genes_combination_numeric[,1:766], 2, ratio_genes_fun, n=NN) # Ratio

genes_name_ratio<-ratio_genes_name_fun(genes_name_combination,n=NN) # Ratio
genes_ratio_all<-as.data.frame(genes_ratio_numeric)
genes_ratio_all$LAG_diff<-LAG_diff_fun(genes_combination_numeric[,767],n=NN)
row.names(genes_ratio_all)<-genes_name_ratio

genes_ratio_all<-genes_ratio_all[which(genes_ratio_all$LAG_diff>=6), ]

genes_ratio_all<-genes_ratio_all[which((24-genes_ratio_all$LAG_diff)>=6), ]

# write.csv(genes_ratio_all, 'Gene_ratio_MAB_6-18hr.csv')