



################################################################################
### Alignment workflow for the four human pancreatic islet datasets
################################################################################
rm(list=ls())

library(Seurat)
library(Matrix)

# Read in all four input expression matrices
# celseq.data <- read.table("~/data/pancreas_multi_celseq_expression_matrix.txt.gz")
# celseq2.data <- read.table("~/data/pancreas_multi_celseq2_expression_matrix.txt.gz")
# fluidigmc1.data <- read.table("~/data/pancreas_multi_fluidigmc1_expression_matrix.txt.gz")
# smartseq2.data <- read.table("~/data/pancreas_multi_smartseq2_expression_matrix.txt.gz")

setwd('C:/Users/shijusis/OneDrive - Michigan Medicine/Desktop/Shiju_sisobhan/RNA sequencing/sRNA Seq/Integration_example')

celseq.data <- read.table("Pancreas_CellSeq.txt")
celseq2.data <- read.table("Pancreas_cellSeq2.csv")

fluidigmc1.data <- read.csv("Pancreas_Fluidigm.csv")
rownames(fluidigmc1.data)<-fluidigmc1.data[,1]
fluidigmc1.data<-fluidigmc1.data[-1]


# Convert to sparse matrices for efficiency
celseq.data <- as(as.matrix(celseq.data), "dgCMatrix")
celseq2.data <- as(as.matrix(celseq2.data), "dgCMatrix")
fluidigmc1.data <- as(as.matrix(fluidigmc1.data), "dgCMatrix")

dim(celseq.data)
dim(celseq2.data)
dim(fluidigmc1.data)



celseq <- CreateSeuratObject(counts = celseq.data)
VlnPlot(celseq, "nFeature_RNA")
celseq <- subset(celseq, subset = nFeature_RNA > 1750)
VlnPlot(celseq, "nFeature_RNA")
celseq <- NormalizeData(celseq, normalization.method = "LogNormalize", scale.factor = 10000)
celseq <- FindVariableFeatures(celseq, selection.method = "vst", nfeatures = 2000)
celseq <- ScaleData(celseq)
celseq[["tech"]] <- "celseq"
View(celseq@meta.data)


# CEL-Seq2 https://www.cell.com/molecular-cell/fulltext/S1097-2765(09)00641-8
# In subset, use low.thresholds = 2500.
celseq2 <- CreateSeuratObject(counts = celseq2.data)
VlnPlot(celseq2, "nFeature_RNA")
celseq2 <- subset(celseq2, subset = nFeature_RNA > 2500)
VlnPlot(celseq2, "nFeature_RNA")
celseq2 <- NormalizeData(celseq2, normalization.method = "LogNormalize", scale.factor = 10000)
celseq2 <- FindVariableFeatures(celseq2, selection.method = "vst", nfeatures = 2000)
celseq2 <- ScaleData(celseq2)
celseq2[["tech"]] <- "celseq2"
View(celseq2@meta.data)


# Fluidigm C1
# Omit subset function because cells are already high quality.
fluidigmc1 <- CreateSeuratObject(counts = fluidigmc1.data)
VlnPlot(fluidigmc1, "nFeature_RNA")
fluidigmc1 <- NormalizeData(fluidigmc1, normalization.method = "LogNormalize", scale.factor = 10000)
fluidigmc1 <- FindVariableFeatures(fluidigmc1, selection.method = "vst", nfeatures = 2000)
fluidigmc1 <- ScaleData(fluidigmc1)
fluidigmc1[["tech"]] <- "fluidigmc1"
View(fluidigmc1@meta.data)


# This code sub-samples the data in order to speed up calculations and not use too much memory.
Idents(celseq) <- "tech"
celseq <- subset(celseq, downsample = 500, seed = 1)
View(celseq@meta.data)

Idents(celseq2) <- "tech"
celseq2 <- subset(celseq2, downsample = 500, seed = 1)

Idents(fluidigmc1) <- "tech"
fluidigmc1 <- subset(fluidigmc1, downsample = 500, seed = 1)


# Merge Seurat objects. Original sample identities are stored in gcdata[["tech"]].
# Cell names will now have the format tech_cellID (smartseq2_cell1...)
add.cell.ids <- c("celseq", "celseq2", "fluidigmc1")
gcdata <- merge(x = celseq, y = list(celseq2, fluidigmc1), add.cell.ids = add.cell.ids, merge.data = FALSE)
Idents(gcdata) <- "tech"  # use identity based on sample identity
View(gcdata@meta.data)

# Look at how the number of genes per cell varies across the different technologies.
VlnPlot(gcdata, "nFeature_RNA", group.by = "tech")


# The merged data must be normalized and scaled (but you only need to scale the variable genes). 
# Let us also find the variable genes again this time using all the pancreas data.
gcdata <- NormalizeData(gcdata, normalization.method = "LogNormalize", scale.factor = 10000)
var.genes <- SelectIntegrationFeatures(SplitObject(gcdata, split.by = "tech"),
                                       nfeatures = 2000, verbose = TRUE, fvf.nfeatures = 2000, selection.method = "vst")
