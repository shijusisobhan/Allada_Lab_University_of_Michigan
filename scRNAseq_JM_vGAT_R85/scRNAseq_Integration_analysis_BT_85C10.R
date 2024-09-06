
rm(list=ls())

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)
library(data.table)


# Read the count matrix for all samples
sparse_matrix <- readMM("X:/Sequencing_data/Clark_seq_data_6-12-2024/01.RawData/BT_all/KB_BT_all/counts_unfiltered/cells_x_genes.mtx")
sparse_matrix<-t(sparse_matrix)
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
BT.Data<-as(Count_matrix, "dgCMatrix")


counts_85C10_ZT0 <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11392-JM/10x_analysis_11392-JM/Sample_11392-JM-1/filtered_feature_bc_matrix/")

counts_85C10_ZT12 <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11399-JM/10x_analysis_11399-JM/Sample_11399-JM-1/filtered_feature_bc_matrix/")

dim(BT.Data)
dim(counts_85C10_ZT0)
dim(counts_85C10_ZT12)


#**************************************************************************************************************************

# Create Seurat Objects for Each Sample

# 1. CreateSeuratObject

obj_BT<-CreateSeuratObject(counts = BT.Data, project="BT")
obj_85C10_ZT0<-CreateSeuratObject(counts = counts_85C10_ZT0, project ='85C10.ZT0' )
obj_85C10_ZT12<-CreateSeuratObject(counts = counts_85C10_ZT12, project = '85C10_ZT12')


# Quality control of each samples

# BT samples are High quality sample- no filtering requird

obj_85C10_ZT0 <- subset(obj_85C10_ZT0, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & 
                          nCount_RNA > 6000 & nCount_RNA < 75000)

obj_85C10_ZT12 <- subset(obj_85C10_ZT12, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & 
                           nCount_RNA > 6000 & nCount_RNA < 75000)



# Normalize the data

obj_BT<-NormalizeData(obj_BT)
obj_85C10_ZT0<-NormalizeData(obj_85C10_ZT0)
obj_85C10_ZT12<-NormalizeData(obj_85C10_ZT12)


# Find Variable Features
#Identify differences in the genes expressed across the datasets.

obj_BT  <- FindVariableFeatures(obj_BT)
obj_85C10_ZT0  <- FindVariableFeatures(obj_85C10_ZT0)
obj_85C10_ZT12  <- FindVariableFeatures(obj_85C10_ZT12)


# Select Integration Features 
# Different data set have different  number of genes.
#To integrate datasets, Seurat requires identifying a common set of features (genes).

features<-SelectIntegrationFeatures(object.list = list(obj_BT,obj_85C10_ZT0,obj_85C10_ZT12))


# Scale the data using the common features

obj_BT<-ScaleData(obj_BT, features = features)
obj_85C10_ZT0<-ScaleData(obj_85C10_ZT0, features = features)
obj_85C10_ZT12<-ScaleData(obj_85C10_ZT12, features = features)




dim(obj_BT)
dim(obj_85C10_ZT0)
dim(obj_85C10_ZT12)


# Perform Integration

# find integration anchors which align the datasets
# Here maximum number of cell is 6 in BT. So that st dim and k.score should be less than 6
anchors<-FindIntegrationAnchors(object.list = list(obj_BT,obj_85C10_ZT0,obj_85C10_ZT12),
                                anchor.features = features, dims = 1:5,  k.score = 5) 



# Integrate the dataset the datasets.
# Here maximum number of cell is 6 in BT. So that k.weight should be less than 6
integrated_data <- IntegrateData(anchorset = anchors, k.weight = 5)


# Further analysis
# scale the integrated data
integrated_data <-ScaleData(object=integrated_data )

# PCA
integrated_data<-RunPCA(object=integrated_data)
DimPlot(integrated_data, reduction = "pca")

# Umap reduction
integrated_data<-RunUMAP(object = integrated_data, dims = 1:20)
plot1<-UMAPPlot(integrated_data)


# Cluster the cells
integrated_data <- FindNeighbors(integrated_data, dims = 1:20)
integrated_data <- FindClusters(integrated_data, resolution = 1)
plot2 <- DimPlot(integrated_data, reduction = "umap", label = T)

plot1+plot2

View(integrated_data@meta.data)




# Find Differentially expressed genes (cluster marker identification)

# Running the IntegrateData function creates a new Assay object (by default it is called integrated)
#The uncorrected values -tore in the original Assay object (called RNA by default).
# The default assay of the resulted Seurat object is automatically set to integrated
# he corrected values are no longer very reliable
#For cluster marker identification and visualization, to use the uncorrected expression values

DefaultAssay(integrated_data) <- "RNA"

# Once integrative analysis is complete, you can rejoin the layers - 
# which collapses the individual datasets together and 
# recreates the original counts and data layers. 

seurat.integrated_join<-JoinLayers(integrated_data)
View(seurat.integrated_join@meta.data)




# find all markers

DEG_all<-FindAllMarkers(seurat.integrated_join,
                        logfc.threshold = 1,
                        min.pct = 0.1,
                        only.pos = T # Only up-regulated genes
)  
# Save the DEG results
write.csv(DEG_all,'X:/Sequencing_data/scRNA_seq_JM_8_14-2024/Data_analysis_Shiju/DEG_Between_cluster_BT_85C10_Integrated.csv')



# Retrieve data in an expression matrix RNA counts matrix
count_matrix_integrated <- as.data.frame(seurat.integrated_join[["RNA"]]$counts)
count_matrix_integrated$genes<-rownames(count_matrix_integrated)

MK1<-FeaturePlot(seurat.integrated_join, features = "Pdf", min.cutoff = 'q10')
MK2<-FeaturePlot(seurat.integrated_join, features = "per", min.cutoff = 'q10', label = T)
MK3<-FeaturePlot(seurat.integrated_join, features = "Clk", min.cutoff = 'q10')
MK4<-FeaturePlot(seurat.integrated_join, features = "tim", min.cutoff = 'q10')
MK5<-FeaturePlot(seurat.integrated_join, features = "vri", min.cutoff = 'q10')
MK6<-FeaturePlot(seurat.integrated_join, features = "AstC", min.cutoff = 'q10')
MK7<-FeaturePlot(seurat.integrated_join, features = "EGFP", min.cutoff = 'q10', label = T)
library(gridExtra)
grid.arrange(MK2,MK4,MK3,MK5,MK1,MK6, ncol=3)



# Perform DE analysis within the same cell type across conditions ************************************************

# 1. create a column in the meta.data slot to hold both the cell type and ZT information



seurat.integrated_join@meta.data$cell_condition <- ifelse(grepl("ZT0", seurat.integrated_join@meta.data$orig.ident),
                                                          paste(seurat.integrated_join@meta.data$seurat_clusters, "ZT0",sep = "_"),
                                                          ifelse(grepl("ZT12", seurat.integrated_join@meta.data$orig.ident),
                                                                 paste(seurat.integrated_join@meta.data$seurat_clusters, "ZT12",sep = "_"),
                                                                 paste(seurat.integrated_join@meta.data$seurat_clusters, "ZT0",sep = "_")))  # Optional: BT, ZT0 for any other cases




View(seurat.integrated_join@meta.data)
# set cell_condition as objects identity classes 
Idents(seurat.integrated_join)<-'cell_condition'

# There are 16 clusters 0-15
DEG_ZT0vsZT12_all<-data.frame()
for (i in 0:15 ) {
  cluster_number<-i
  DEG_ZT0vsZT12<- FindMarkers(seurat.integrated_join, ident.1 = paste(cluster_number,"ZT0", sep ="_"), 
                              ident.2 = paste(cluster_number,"ZT12", sep ="_"), verbose = FALSE)
  DEG_ZT0vsZT12$cluster<-cluster_number
  DEG_ZT0vsZT12$Genes<-rownames(DEG_ZT0vsZT12)
  
  DEG_ZT0vsZT12_all<-rbind(DEG_ZT0vsZT12_all,DEG_ZT0vsZT12[which(DEG_ZT0vsZT12$p_val_adj<0.1),])
}


write.csv(DEG_ZT0vsZT12_all,'X:/Sequencing_data/scRNA_seq_JM_8_14-2024/Data_analysis_Shiju/DEG_Between_condition_BT_85C10_ZT0vsZT12.csv')
