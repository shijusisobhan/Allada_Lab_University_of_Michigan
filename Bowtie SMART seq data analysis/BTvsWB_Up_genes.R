
rm(list=ls())

setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis')

DEG_all<-read.csv("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BT_allvsWB_ZT0.csv")
DEG_cell1<-read.csv("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BT_cell1vsWB_ZT0.csv")
DEG_cell2<-read.csv("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BT_cell2vsWB_ZT0.csv")
DEG_cell3<-read.csv("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BT_cell3vsWB_ZT0.csv")
DEG_cell4<-read.csv("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BT_cell4vsWB_ZT0.csv")
DEG_cell5<-read.csv("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BT_cell5vsWB_ZT0.csv")
DEG_cell6<-read.csv("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BT_cell6vsWB_ZT0.csv")


columns_list<-list(
  
  gene_all<-DEG_all[which(DEG_all$log2FoldChange>0.6 & DEG_all$padj<0.05), 'ext_gene'],
  gene_cell1<-DEG_cell1[which(DEG_cell1$log2FoldChange>0.6 & DEG_cell1$padj<0.05), 'X'],
  gene_cell2<-DEG_cell2[which(DEG_cell2$log2FoldChange>0.6 & DEG_cell2$padj<0.05), 'X'],
  gene_cell3<-DEG_cell3[which(DEG_cell3$log2FoldChange>0.6 & DEG_cell3$padj<0.05), 'X'],
  gene_cell4<-DEG_cell4[which(DEG_cell4$log2FoldChange>0.6 & DEG_cell4$padj<0.05), 'X'],
  gene_cell5<-DEG_cell5[which(DEG_cell5$log2FoldChange>0.6 & DEG_cell5$padj<0.05), 'X'],
  gene_cell6<-DEG_cell6[which(DEG_cell6$log2FoldChange>0.6 & DEG_cell6$padj<0.05), 'X']
)


# Find the maximum length of the columns
max_length <- max(sapply(columns_list, length))

# Pad shorter columns with NA to match the maximum length
columns_list <- lapply(columns_list, function(col) {
  c(col, rep(NA, max_length - length(col)))
})

# Combine columns into a data frame
df <- data.frame(columns_list)

colnames(df)<- c('BT_all','BT1', 'BT2','BT3','BT4','BT5','BT6')

# Write the data frame to a CSV file
write.csv(df, "BTvsWB_Upregulated_genes.csv", row.names = FALSE)


