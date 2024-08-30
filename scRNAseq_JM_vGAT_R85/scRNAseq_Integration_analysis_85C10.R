rm(list=ls())

library(Seurat)
library(dplyr)

counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11392-JM/10x_analysis_11392-JM/Sample_11392-JM-1/filtered_feature_bc_matrix/")
DBD.85C10.VT.AD_ZT0 <- CreateSeuratObject(counts, project="DBD.85C10.VT.AD_ZT0")


counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11399-JM/10x_analysis_11399-JM/Sample_11399-JM-1/filtered_feature_bc_matrix/")
DBD.85C10.VT.AD_ZT12 <- CreateSeuratObject(counts, project="DBD.85C10.VT.AD_ZT12")


# Merge data set (Not Intgrete) and do quality control (ls() function give you the list of objects)

merged_seurat<-merge(DBD.85C10.VT.AD_ZT0,y=DBD.85C10.VT.AD_ZT12, add.cell.ids=ls()[2:3], project='JM')

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
top_features <- head(VariableFeatures(merged_seurat_filtered ), 10)
plot1 <- VariableFeaturePlot(merged_seurat_filtered )
# LabelPoints(plot = plot1, points = top_features, repel = TRUE)


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

plot1 + plot3
plot2+plot4



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

plot5 <- TSNEPlot(seurat.integrated)
plot6 <- UMAPPlot(seurat.integrated)
plot5 + plot6



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
write.csv(DEG_all,'X:/Sequencing_data/scRNA_seq_JM_8_14-2024/Data_analysis_Shiju/DEG_Between_cluster_85C10_all.csv')

# Let’s take a quick glance at the markers.
# find out number of up-regulated genes in each cluster compared to other clusters
table(DEG_all$cluster)

# find top 3 genes in each clusters based on log2FC

top3_markers <- as.data.frame(DEG_all %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
top3_markers

MK1<-FeaturePlot(seurat.integrated_join, features = "AstC", min.cutoff = 'q10')
MK2<-FeaturePlot(seurat.integrated_join, features = "twit", min.cutoff = 'q10')
MK3<-FeaturePlot(seurat.integrated_join, features = "C15", min.cutoff = 'q10')
MK4<-FeaturePlot(seurat.integrated_join, features = "EGFP", min.cutoff = 'q10')
library(gridExtra)
grid.arrange(MK4,MK1,MK2,MK3, ncol=2)

#***********************************************************************************
# Find number of cells and genes

num_cells <- ncol(seurat.integrated_join)
num_features <- nrow(seurat.integrated_join)

# Retrieve data in an expression matrix RNA counts matrix
count_matrix_integrated <- as.data.frame(seurat.integrated_join[["RNA"]]$counts)
count_matrix_integrated$genes<-rownames(count_matrix_integrated)

MK1<-FeaturePlot(seurat.integrated_join, features = "Pdf", min.cutoff = 'q10')
MK2<-FeaturePlot(seurat.integrated_join, features = "per", min.cutoff = 'q10')
MK3<-FeaturePlot(seurat.integrated_join, features = "Clk", min.cutoff = 'q10')
MK4<-FeaturePlot(seurat.integrated_join, features = "tim", min.cutoff = 'q10')
MK5<-FeaturePlot(seurat.integrated_join, features = "vri", min.cutoff = 'q10')
MK6<-FeaturePlot(seurat.integrated_join, features = "AstC", min.cutoff = 'q10')
library(gridExtra)
grid.arrange(MK2,MK4,MK3,MK5,MK1,MK6, ncol=3)

# Perform DE analysis within the same cell type across conditions (ZT0 vs ZT12) ************************************************

# 1. create a column in the meta.data slot to hold both the cell type and ZT information

seurat.integrated_join@meta.data$cell_condition <- paste(seurat.integrated_join@meta.data$seurat_clusters, seurat.integrated_join@meta.data$orig.ident, sep = "_")
View(seurat.integrated_join@meta.data)
Idents(seurat.integrated_join)<-'cell_condition'

# There are 22 clusters 0-21
DEG_ZT0vsZT12_all<-data.frame()
for (i in 0:21 ) {
  cluster_number<-i
  DEG_ZT0vsZT12<- FindMarkers(seurat.integrated_join, ident.1 = paste(cluster_number,"DBD.85C10.VT.AD_ZT0", sep ="_"), 
                              ident.2 = paste(cluster_number,"DBD.85C10.VT.AD_ZT12", sep ="_"), verbose = FALSE)
  DEG_ZT0vsZT12$cluster<-cluster_number
  DEG_ZT0vsZT12$Genes<-rownames(DEG_ZT0vsZT12)
  
  DEG_ZT0vsZT12_all<-rbind(DEG_ZT0vsZT12_all,DEG_ZT0vsZT12[which(DEG_ZT0vsZT12$p_val_adj<0.1),])
}


write.csv(DEG_ZT0vsZT12_all,'X:/Sequencing_data/scRNA_seq_JM_8_14-2024/Data_analysis_Shiju/DEG_between_condition_85C10_ZT0vsZT12.csv')

# ********************************************************************************************************************

# *****************pseudobulk analysis ********************************

# Above tests treat each cell as an independent replicate and ignore 
# inherent correlations between cells originating from the same sample.
# So, large number of false positive 


# pseudobulking steps
# sum together gene counts of all the cells from the same sample for each cells type (cluster)
# This results in one gene expression profile per sample and cell type (cluster)
# perform DE analysis using DESeq2 on the sample level. 
# This treats the samples, rather than the individual cells, as independent observations.

# Pseudo bulk analysis is usually performed before Integrating the the data

View(merged_seurat_filtered@meta.data)

# aggregate the counts in a samples (orig.ident) with similar cell-cluster (seurat_clusters)

pseudo_bulk_obj<-AggregateExpression(merged_seurat_filtered,assays = "RNA", 
                                     return.seurat = F, group.by = c("seurat_clusters","orig.ident"),
                                     slot='counts')


# Now view the aggregated count matrix 

Pseudo_bulk_count_matrix <- as.matrix(pseudo_bulk_obj$RNA)
Pseudo_bulk_count_matrix[1:5,1:5] 
# within cell cluster-0 (g0), sample ZT0 (DBD.85C10.VT.AD-ZT0) 
# has 135 counts associated with gene TyrR.

# Now rows- genes, cloumn- samples
# we need to split the data according to cluster (cell type) (for DEG analysis)
# To do that first transpose the matrix such that, rows-samples, columns- genes

Pseudo_bulk_count_matrix_t<-as.data.frame(t(Pseudo_bulk_count_matrix))

cluster_name <- sub("_.*", "", rownames(Pseudo_bulk_count_matrix_t)) #Replaces only the first occurrence
print(cluster_name) # 17 clusters before integration

# split according to cluster name
Pseudobulk_split<-split.data.frame(Pseudo_bulk_count_matrix_t,
                                   f=factor(cluster_name))

Pseudobulk_split$g0[1:2,1:5]

# Now we have to do two more steps
# 1. Exclude cluster name from rows (only need sample name (for bulk rna seq))
# 2. transpose back the count matrix (rows-genes, column-samples)

Pseudobulk_split_modified<-lapply(Pseudobulk_split, function(x){
  rownames(x) <- sub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

Pseudobulk_split_modified$g0[1:2, 1:2]

# Now prepare data for DEseq2 analysis
# Start with data in cluster 0
# 1. get the count matrix 
# 2. Create the sample condition table (metadata)

counts_g0<-Pseudobulk_split_modified$g0

# DEseq2 don't work here because there is only two samples.


