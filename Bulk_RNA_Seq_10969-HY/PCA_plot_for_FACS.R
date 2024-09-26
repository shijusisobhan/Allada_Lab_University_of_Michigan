rm(list=ls())

C_dir<- "X:/Sequencing_data/HY_10969/Data_analysis_shiju"
setwd(C_dir)

library(DESeq2)
rawCounts <- read.csv("Estimated_count_FACS.csv")

# Read in the sample mappings
sampleData <- read.csv('sample_condition_FACS.csv')

geneID <- rawCounts$ext_gene
rawCounts<-rawCounts[-1]
rownames(rawCounts) <- geneID
rawCounts <- as.matrix(rawCounts)
rawCounts<-round(rawCounts, digits=0)

colnames(rawCounts)<-sampleData$sample
rawCounts<-na.omit(rawCounts)

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData) <- sampleData$sample

sampleData <- sampleData[-1]

Condition<-unique(sampleData$Condition)

sampleData$Condition <- factor(sampleData$Condition, 
                            levels=Condition)

# Create the DEseq2DataSet object
deseq2Data <- DESeqDataSetFromMatrix(countData=rawCounts, colData=sampleData, design= ~ Condition)

# Perform pre-filtering of the data
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 5, ]


deseq2Data <- DESeq(deseq2Data)

vsd <- vst(deseq2Data, blind=FALSE)

plotPCA(vsd, intgroup = "Condition")