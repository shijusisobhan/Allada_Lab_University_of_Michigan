rm(list=ls())

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)
library(data.table)

# Read the count matrix for bowtie nuerons (ZT0) single cell ena seq data
sparse_matrix <- readMM("X:/Sequencing_data/Clark_seq_data_6-12-2024/01.RawData/BT_all/KB_BT_all/counts_unfiltered/cells_x_genes.mtx")
sparse_matrix<-t(sparse_matrix)

# Read the genes and barcodes files
genes <- fread("X:/Sequencing_data/Clark_seq_data_6-12-2024/01.RawData/BT_all/KB_BT_all/counts_unfiltered/cells_x_genes.genes.names.txt", header = FALSE)
barcodes <- fread("X:/Sequencing_data/Clark_seq_data_6-12-2024/01.RawData/BT_all/KB_BT_all/counts_unfiltered/cells_x_genes.barcodes.txt", header = FALSE)

# Convert to vectors
gene_names <- genes$V1
barcode_names <- barcodes$V1

# Create a vector of meaningful names
meaningful_names <- paste0("Cell", seq_along(barcode_names))

# Assign row names (genes) and column names (barcodes) to the matrix
rownames(sparse_matrix) <- gene_names
colnames(sparse_matrix)<-meaningful_names
Count_matrix <- as.matrix(sparse_matrix)
Count_matrix <- Count_matrix[rowSums(Count_matrix) > 0, ]


## ********************************************************************************************************    

# Normalize Bulk RNA-seq Data (DESeq2):

bulk_data <- read.csv("C:\\Users\\shijusis\\OneDrive - Michigan Medicine\\Desktop\\Shiju_sisobhan\\RNA sequencing\\Drosophila\\Brain data\\Clark data brain\\Strica_C\\Circadian_profiling\\Estimated_counts_CR_profiling.csv")
rawCounts<-bulk_data[,1:4]

rawCounts_bulk<-rawCounts[-1]
rownames(rawCounts_bulk)<-rawCounts[,1]
rawCounts_bulk <- rawCounts_bulk[rowSums(rawCounts_bulk) > 0, ]



#*********************************************************************************************
# combined individual cell data and all the bulk data for further analysis

combined_data<-merge(Count_matrix,rawCounts_bulk,by = "row.names", all = F)

#*********************************************************************************************  
#*
library(DESeq2)

# Now do the DGE analysis using Deseq2 with individual cells
# for all cell

# geneID_BT_WB<-combined_data[1]
# rawCounts_BT_WB<-combined_data[-1]
# sampleData_BT_WB<-data.frame(sample=(colnames(rawCounts_BT_WB)), condition=c(rep("BT",6), rep("WB",3)))

# *************************************************************************************************************

# for one cell

rawCounts_BT_WB<-combined_data[-c(2:6)] # for one cell
# rawCounts_BT_WB<-combined_data[-c(2:6)] # for one cell

## Exclude zero conts from each individual BT cells###################
rawCounts_BT_WB<-rawCounts_BT_WB[rawCounts_BT_WB[,2] !=0,] # for one cell
geneID_BT_WB<-rawCounts_BT_WB[1]
rawCounts_BT_WB<-rawCounts_BT_WB[-1]


sampleData_BT_WB<-data.frame(sample=(colnames(rawCounts_BT_WB)), condition=c(rep("BT",1), rep("WB",3)))

# *****************************************************************************************
rownames(rawCounts_BT_WB) <- geneID_BT_WB$Row.names
rawCounts_BT_WB <- as.matrix(rawCounts_BT_WB)
rawCounts_BT_WB<-round(rawCounts_BT_WB, digits=0)
colnames(rawCounts_BT_WB)<-sampleData_BT_WB$sample
rawCounts_BT_WB<-na.omit(rawCounts_BT_WB)

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData_BT_WB) <- sampleData_BT_WB$sample

sampleData_BT_WB <- sampleData_BT_WB[-1]

sampleData_BT_WB$condition <- factor(sampleData_BT_WB$condition, 
                                     levels=c("BT","WB"))
# Create the DEseq2DataSet object
deseq2Data_BT_WB <- DESeqDataSetFromMatrix(countData=rawCounts_BT_WB, colData=sampleData_BT_WB, design= ~condition)
deseq2Data_BT_WB <- DESeq(deseq2Data_BT_WB)
vsd <- vst(deseq2Data_BT_WB, blind=FALSE)

plotPCA(vsd)


#deseq2Results_BT_WB <- results(deseq2Data_BT_WB)

deseq2Results_BT_WB <- results(deseq2Data_BT_WB, contrast=c('condition','BT', 'WB'))

# contrast = c('factorName','numeratorLevel','denominatorLevel')


deseq2Results_BT_WB<-na.omit(deseq2Results_BT_WB)

resOrdered <- deseq2Results_BT_WB[order(deseq2Results_BT_WB$pvalue),]

# Normalized_sigle_bulk <- counts(deseq2Data_BT_WB, normalized = TRUE)
# write.csv(as.data.frame(Normalized_sigle_bulk),
#           file="X:/Sequencing_data/Clark_seq_data_6-12-2024/01.RawData/BT_all/Normalized_data_sigle_bulk.csv")

test_table<-as.data.frame(resOrdered)
test_table$ext_gene<-row.names(test_table)

test_table$neg_log10_qval<- -log10(test_table$padj)

#********************************************************************************************
# To exclude neg_log10_qval =inf
finite_values<-test_table$neg_log10_qval[is.finite(test_table$neg_log10_qval)]

max_value<-max(finite_values) #Set a Cap on the Maximum Value
test_table$neg_log10_qval<-pmin(test_table$neg_log10_qval,max_value)
#****************************************************************************** *************


sig_level<-0.05

test_table$diffexpressed <- "Not sig"
# if log2Foldchange > 0 and qvalue < 0.1, set as "UP" 
test_table$diffexpressed[test_table$log2FoldChange > 1 & test_table$padj < sig_level] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
test_table$diffexpressed[test_table$log2FoldChange < -1 & test_table$padj < sig_level] <- "DOWN"

N_significant<-length(test_table$diffexpressed[test_table$diffexpressed !="Not sig"])
N_UP<-length(test_table$diffexpressed[test_table$diffexpressed =="UP"])
N_DOWN<-length(test_table$diffexpressed[test_table$diffexpressed =="DOWN"])

N_significant
N_UP
N_DOWN

UP_genes<-test_table[test_table$diffexpressed=="UP",]
Down_genes<-test_table[test_table$diffexpressed=="DOWN",]
top_genes<- UP_genes[1:3,]

library('ggplot2')
library(ggrepel)
ggplot(test_table) + geom_point(aes(x = log2FoldChange, y = neg_log10_qval, col=diffexpressed))+
  geom_vline(xintercept=0, col="black",  linetype="dashed")+
  scale_color_manual(values=c("blue", "black", "red"))+
  labs(x = "log2(FC)", y="-log10(q)", colour="DEG") +
  geom_text_repel(data = top_genes, 
                  aes(x = log2FoldChange, y = neg_log10_qval, label = ext_gene), color = "black",
                  min.segment.length = unit(0, 'lines'), nudge_y = 20)

write.csv(test_table[which(test_table$padj<0.05 & abs(test_table$log2FoldChange)>1),], 
          file="C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BT6vsWB_ZT0.csv")



