rm(list=ls())

library(Seurat)
library(dplyr)

# # Create a Seurat object
counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11358-JM/10x_analysis_11358-JM/Sample_11358-JM-1/filtered_feature_bc_matrix/")
vGAT_ZT0_1 <- CreateSeuratObject(counts, project="vGAT_ZT0_1")

counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11358-JM/10x_analysis_11358-JM/Sample_11358-JM-2/filtered_feature_bc_matrix/")
vGAT_ZT0_2 <- CreateSeuratObject(counts, project="vGAT_ZT0_2")

counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11363-JM/10x_analysis_11363-JM/Sample_11363-JM-1/filtered_feature_bc_matrix/")
vGAT_ZT12_1 <- CreateSeuratObject(counts, project="vGAT_ZT12_1")

counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11363-JM/10x_analysis_11363-JM/Sample_11363-JM-2/filtered_feature_bc_matrix/")
vGAT_ZT12_2 <- CreateSeuratObject(counts, project="vGAT_ZT12_2")


# Merge data set (Not Intgrete) and do quality control (ls() function give you the list of objects)

merged_seurat<-merge(vGAT_ZT0_1,y=c(vGAT_ZT0_2,vGAT_ZT12_1,vGAT_ZT12_2), add.cell.ids=ls()[2:5], project='JM_vGAT')

View(merged_seurat@meta.data)

merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT[-\\.]")

VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Exclude the cell with
#1. fewer than 1000 or more than 6000 genes
# 2. fewer than 6000 or more than 75000 total UMI
# 3. mt percentage >5%

merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & 
                                   nCount_RNA > 6000 & nCount_RNA < 75000 &
                                   percent.mt < 5)


# Normalization
merged_seurat_filtered <- NormalizeData(merged_seurat_filtered)


merged_seurat_filtered  <- FindVariableFeatures(merged_seurat_filtered , nfeatures = 3000)
top_features <- head(VariableFeatures(merged_seurat_filtered ), 3)
plot1 <- VariableFeaturePlot(merged_seurat_filtered )
#LabelPoints(plot = plot1, points = top_features, repel = TRUE)


# scale data

merged_seurat_filtered <- ScaleData(object=merged_seurat_filtered)

# PCA

merged_seurat_filtered <- RunPCA(merged_seurat_filtered, npcs = 50)
ElbowPlot(merged_seurat_filtered, ndims = ncol(Embeddings(merged_seurat_filtered, "pca")))
# Plot the cells in the 2D PCA projection
DimPlot(merged_seurat_filtered, reduction = "pca")


# Non-linear dimension reduction for visualization
merged_seurat_filtered <- RunTSNE(merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:20)



# Cluster the cells
merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(merged_seurat_filtered, resolution = 1)

plot1 <- DimPlot(merged_seurat_filtered, reduction = "tsne", label = F)
plot2 <- DimPlot(merged_seurat_filtered, reduction = "umap", label = F)


# by grouping examining the batch effect
plot3 <- DimPlot(merged_seurat_filtered, reduction = "tsne", label = F, group.by = 'orig.ident')
plot4 <- DimPlot(merged_seurat_filtered, reduction = "umap", label = F, group.by ='orig.ident')

plot3 + plot1
plot4+plot2



# Now integrate data

# aim to integrate data from the two conditions, so that 
# cells from the same cell type/subpopulation will cluster together.

obj.list<-SplitObject(merged_seurat_filtered, split.by = 'orig.ident')

for (i in 1:length(obj.list)) {
  
  obj.list[[i]]<-NormalizeData(object=obj.list[[i]])
  obj.list[[i]]<-FindVariableFeatures(object=obj.list[[i]])
  
}

# select integration features

features<-SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors
anchors<-FindIntegrationAnchors(object.list = obj.list,
                                anchor.features = features)

# integrate data
seurat.integrated<-IntegrateData(anchorset = anchors)

# scale the data

seurat.integrated<-ScaleData(object=seurat.integrated)

# Run PCA

seurat.integrated<-RunPCA(object=seurat.integrated)
DimPlot(seurat.integrated, reduction = "pca")

#Run umap

seurat.integrated<-RunUMAP(object = seurat.integrated, dims = 1:20)
seurat.integrated <- RunTSNE(seurat.integrated, dims = 1:20)


# Cluster the cells
seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:20)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 1)

plot1 <- DimPlot(seurat.integrated, reduction = "tsne", label = F)
plot2 <- DimPlot(seurat.integrated, reduction = "umap", label = F)



# by grouping examining the batch effect
plot3 <- DimPlot(seurat.integrated, reduction = "tsne", label = F, group.by = 'orig.ident')
plot4 <- DimPlot(seurat.integrated, reduction = "umap", label = F, group.by ='orig.ident')

plot3 + plot1
plot4+plot2

View(seurat.integrated@meta.data)

# Find Differentially expressed genes (cluster marker identification)

# Running the IntegrateData function creates a new Assay object (by default it is called integrated)
#The uncorrected values -tore in the original Assay object (called RNA by default).
# The default assay of the resulted Seurat object is automatically set to integrated
# he corrected values are no longer very reliable
#For cluster marker identification and visualization, to use the uncorrected expression values

DefaultAssay(seurat.integrated) <- "RNA"

# Once integrative analysis is complete, you can rejoin the layers - 
# which collapses the individual datasets together and 
# recreates the original counts and data layers. 

seurat.integrated_join<-JoinLayers(seurat.integrated)
View(seurat.integrated_join@meta.data)




# find all markers

DEG_all<-FindAllMarkers(seurat.integrated_join,
                        logfc.threshold = 1,
                        min.pct = 0.1,
                        only.pos = T # Only up-regulated genes
)  
# Save the DEG results
write.csv(DEG_all,'X:/Sequencing_data/scRNA_seq_JM_8_14-2024/Data_analysis_Shiju/DEG_Between_cluster_vGAT_ZT0_ZT12.csv')

# Let’s take a quick glance at the markers.
# find out number of up-regulated genes in each cluster compared to other clusters
table(DEG_all$cluster)

# find top 3 genes in each clusters based on log2FC

top3_markers <- as.data.frame(DEG_all %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers

MK1<-FeaturePlot(seurat.integrated_join, features = "EGFP", min.cutoff = 'q10')
MK2<-FeaturePlot(seurat.integrated_join, features = "Imp", min.cutoff = 'q10')
MK3<-FeaturePlot(seurat.integrated_join, features = "AstA", min.cutoff = 'q10')
MK4<-FeaturePlot(seurat.integrated_join, features = "AstC", min.cutoff = 'q10')
library(gridExtra)
grid.arrange(MK1,MK2,MK3,MK4, ncol=2)

FeaturePlot(seurat.integrated_join, features = "EGFP", min.cutoff = 'q10', label = T)

# Perform DE analysis within the same cell type across conditions ************************************************

# 1. create a column in the meta.data slot to hold both the cell type and ZT information

#seurat.integrated_join@meta.data$cell_condition <- paste(seurat.integrated_join@meta.data$seurat_clusters, seurat.integrated_join@meta.data$orig.ident, sep = "_")



seurat.integrated_join@meta.data$cell_condition <- ifelse(grepl("ZT0_", seurat.integrated_join@meta.data$orig.ident),
                                                          paste(seurat.integrated_join@meta.data$seurat_clusters, "vGAT_ZT0",sep = "_"),
                                                          ifelse(grepl("ZT12", seurat.integrated_join@meta.data$orig.ident),
                                                                 paste(seurat.integrated_join@meta.data$seurat_clusters, "vGAT_ZT12",sep="_"),
                                                                 ifelse(grepl("ZT0SD", seurat.integrated_join@meta.data$orig.ident),
                                                                        paste(seurat.integrated_join@meta.data$seurat_clusters, "vGAT_ZT0SD",sep="_"),
                                                                        NA)))  # Optional: NA for any other cases



View(seurat.integrated_join@meta.data)
# set cell_condition as objects identity classes 
Idents(seurat.integrated_join)<-'cell_condition'


# There are 41 clusters 0-20
DEG_ZT0vsZT12_all<-data.frame()
for (i in 0:40 ) {
  cluster_number<-i
  DEG_ZT0vsZT12<- FindMarkers(seurat.integrated_join, ident.1 = paste(cluster_number,"vGAT_ZT0", sep ="_"), 
                              ident.2 = paste(cluster_number,"vGAT_ZT12", sep ="_"), verbose = FALSE)
  DEG_ZT0vsZT12$cluster<-cluster_number
  DEG_ZT0vsZT12$Genes<-rownames(DEG_ZT0vsZT12)
  
  DEG_ZT0vsZT12_all<-rbind(DEG_ZT0vsZT12_all,DEG_ZT0vsZT12[which(DEG_ZT0vsZT12$p_val_adj<0.1),])
}


write.csv(DEG_ZT0vsZT12_all,'X:/Sequencing_data/scRNA_seq_JM_8_14-2024/Data_analysis_Shiju/DEG_Between_condition_vGAT_ZT0vsZT12.csv')

# ********************************************************************************************************************
# ************ Perform DE analysis after pseudobulking****************************************

# Above tests treat each cell as an independent replicate and ignore 
# inherent correlations between cells originating from the same sample.
# So, large number of false positive 


  # pseudobulking steps
# sum together gene counts of all the cells from the same sample for each cellstype
# This results in one gene expression profile per sample and cell type
# perform DE analysis using DESeq2 on the sample level. 
# This treats the samples, rather than the individual cells, as independent observations.

