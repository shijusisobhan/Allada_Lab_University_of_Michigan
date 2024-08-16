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
#bulk_data <- read.csv("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/Drosophila/Fat body data/Whole_Brain_Vs_Fatbody/Estimated_counts_MB247.csv")

#bulk_data <- read.csv("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/Drosophila/Fat body data/Whole_Brain_Vs_Fatbody/Estimated_counts_R85_ZT0vsZT12.csv")

#bulk_data <- read.csv("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/Drosophila/Fat body data/Whole_Brain_Vs_Fatbody/Estimated_counts_ME.csv")

bulk_data <- read.csv("C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/Drosophila/Fat body data/Whole_Brain_Vs_Fatbody/Estimated_counts_vGAT.csv")


#rawCounts<-bulk_data[,c(1,8,9,10)] # for R5_MC
#rawCounts<-bulk_data[,c(1,2,8,9)] # for MB247 ZT0 +ZT12
#rawCounts<-bulk_data[,c(1:4)] # for fACS R85 ZT0
#rawCounts<-bulk_data[,c(1:7)] # for fACS R85 ZT0+ZT12
rawCounts<-bulk_data[,c(1,5:10)] # for vGAT ZT0

rawCounts_bulk<-rawCounts[-1]
rownames(rawCounts_bulk)<-rawCounts[,1]
#rawCounts_bulk <- rawCounts_bulk[rowSums(rawCounts_bulk) > 0, ]
rawCounts_bulk <- subset(rawCounts_bulk,rawCounts_bulk>0)


#*********************************************************************************************
# combined individual cell data and all the bulk data for further analysis

combined_data<-merge(Count_matrix,rawCounts_bulk,by = "row.names", all = F)

library(DESeq2)
#*********************************************************************************************  
# Now do the DGE analysis using Deseq2 with individual cells
# for all cell

 geneID_Combined<-combined_data[1]
 rawCounts_Combined<-combined_data[-1]
 
 second_variable<-"vGAT_ZT0+ZT12" # provide second variable name (eg: "MB247")
 
 sampleData_Combined<-data.frame(sample=(colnames(rawCounts_Combined)), condition=c(rep("BT",6), rep(second_variable,6)))

 #sampleData_Combined<-data.frame(sample=(colnames(rawCounts_Combined)), condition=c(rep("BT",6), rep("R5_MC",3)))

 #*******************************************************************************************************
 # for one cell
  #rawCounts_Combined<-combined_data[-c(2:6)] # for one cell
  #rawCounts_Combined<-combined_data[-c(2:6)] # for one cell
 
#  # Exclude zero conts from each individual BT cells
  # rawCounts_Combined<-rawCounts_Combined[rawCounts_Combined[,2] !=0,] # for one cell
  # geneID_Combined<-rawCounts_Combined[1]
  # rawCounts_Combined<-rawCounts_Combined[-1]
  # sampleData_Combined<-data.frame(sample=(colnames(rawCounts_Combined)), condition=c(rep("BT",1), rep(second_variable,3)))
 

 #*************************************************************************************************************

rownames(rawCounts_Combined) <- geneID_Combined$Row.names
rawCounts_Combined <- as.matrix(rawCounts_Combined)
rawCounts_Combined<-round(rawCounts_Combined, digits=0)
colnames(rawCounts_Combined)<-sampleData_Combined$sample
rawCounts_Combined<-na.omit(rawCounts_Combined)

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData_Combined) <- sampleData_Combined$sample

sampleData_Combined <- sampleData_Combined[-1]

sampleData_Combined$condition <- factor(sampleData_Combined$condition, 
                                     levels=c("BT",second_variable))
# Create the DEseq2DataSet object
deseq2Data_Combined <- DESeqDataSetFromMatrix(countData=rawCounts_Combined, colData=sampleData_Combined, design= ~condition)
deseq2Data_Combined <- DESeq(deseq2Data_Combined)
vsd <- vst(deseq2Data_Combined, blind=FALSE)

plotPCA(vsd)


deseq2Results_Combined <- results(deseq2Data_Combined, contrast=c('condition','BT', second_variable)) 
# contrast = c('factorName','numeratorLevel','denominatorLevel')


deseq2Results_Combined<-na.omit(deseq2Results_Combined)

resOrdered <- deseq2Results_Combined[order(deseq2Results_Combined$pvalue),]


test_table<-as.data.frame(resOrdered)
test_table$ext_gene<-row.names(test_table)

test_table$neg_log10_qval<- -log10(test_table$padj)

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
          file="C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Bowtai_SMARTseq_data_analysis/DEG_BTvsvGAT_ZT0_ZT12.csv")


