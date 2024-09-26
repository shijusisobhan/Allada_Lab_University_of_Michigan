
# collect estimate count and TPM for HY samples
rm(list=ls())
library("tximport")

dir_1<-('X:/Sequencing_data/HY_10969')
samples <- read.table(file.path(dir_1,"sample.csv"), header=TRUE)

dir_2<-('X:/Sequencing_data/HY_10969/fastqs_10969-HY/Kallisto_out_SingleEnd')
files<-file.path(dir_2,samples$sample,'abundance.tsv')
names(files)<-paste0("HY_",samples$sample)

dir_3<-('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/Drosophila/Brain data/Clark data brain/Strica_C/SD/NonSDvs3hrSD/DEseq2_test')
txt2gene <- read.csv(file.path(dir_3,"tr2gnName.csv"), header=TRUE)

txi<-tximport(files, type='kallisto', tx2gene=txt2gene)

write.csv(txi[["counts"]],"X:/Sequencing_data/HY_10969/Data_analysis_shiju/Estimated_counts_HY.csv")
write.csv(txi[["abundance"]],"X:/Sequencing_data/HY_10969/Data_analysis_shiju/TPM_HY.csv")

#**********************************************************************************************************

# Now combine estimate count for HY, MB247, DN1, LNv, R5 samples

library(dplyr)
library(readr)

C_dir<- "X:/Sequencing_data/HY_10969/Data_analysis_shiju"
setwd(C_dir)
est_files<-list.files(C_dir,pattern = "^Estimated", full.names = T,ignore.case = T )

est_all<-read.csv(est_files[1])

colnames(est_all)[1]<-'ext_gene'

for (i in 2:length(est_files)) {
  
  D1<-read.csv(est_files[i])
  colnames(D1)[1]<-'ext_gene' # first column name should be 'ext_gene'
  est_all<-merge(est_all, D1,  by.x=1 , by.y="ext_gene" )
  
}

samples_names<-colnames(est_all)

write.csv(samples_names, 'sample_condition_FACS.csv')
write.csv(est_all, 'Estimated_count_FACS.csv')
