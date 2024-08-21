rm(list=ls())

library(Seurat)

# Create a Seurat object
counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11358-JM/10x_analysis_11358-JM/Sample_11358-JM-1/filtered_feature_bc_matrix/")
vGAT_ZT0_1 <- CreateSeuratObject(counts, project="vGAT_ZT0_1")

counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11358-JM/10x_analysis_11358-JM/Sample_11358-JM-2/filtered_feature_bc_matrix/")
vGAT_ZT0_2 <- CreateSeuratObject(counts, project="vGAT_ZT0_2")

counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11363-JM/10x_analysis_11363-JM/Sample_11363-JM-1/filtered_feature_bc_matrix/")
vGAT_ZT12_1 <- CreateSeuratObject(counts, project="vGAT_ZT12_1")

counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11363-JM/10x_analysis_11363-JM/Sample_11363-JM-2/filtered_feature_bc_matrix/")
vGAT_ZT12_2 <- CreateSeuratObject(counts, project="vGAT_ZT12_2")

counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11373-JM/10x_analysis_11373-JM/Sample_11373-JM-1/filtered_feature_bc_matrix/")
vGAT_ZT0SD_1 <- CreateSeuratObject(counts, project="vGAT_ZT0SD_1")

counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11373-JM/10x_analysis_11373-JM/Sample_11373-JM-2/filtered_feature_bc_matrix/")
vGAT_ZT0SD_2 <- CreateSeuratObject(counts, project="vGAT_ZT0SD_2")

counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11392-JM/10x_analysis_11392-JM/Sample_11392-JM-1/filtered_feature_bc_matrix/")
DBD_85C10_VT_AD_ZT0_1 <- CreateSeuratObject(counts, project="DBD_85C10_VT_AD_ZT0_1")


counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11399-JM/10x_analysis_11399-JM/Sample_11399-JM-1/filtered_feature_bc_matrix/")
DBD_85C10_VT_AD_ZT12_1 <- CreateSeuratObject(counts, project="DBD_85C10_VT_AD_ZT12_1")


# Merge data set (Not Intgrete) and do quality control (ls() function give you the list of objects)

merged_seurat<-merge(DBD_85C10_VT_AD_ZT0_1,y=c(DBD_85C10_VT_AD_ZT12_1,vGAT_ZT0_1,vGAT_ZT0_2,vGAT_ZT0SD_1,vGAT_ZT0SD_2,
                                              vGAT_ZT12_1,vGAT_ZT12_2), add.cell.ids=ls()[2:9], project='JM')

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
LabelPoints(plot = plot1, points = top_features, repel = TRUE)


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

plot1 <- TSNEPlot(merged_seurat_filtered)
plot2 <- UMAPPlot(merged_seurat_filtered)
plot1 + plot2


# Cluster the cells
merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(merged_seurat_filtered, resolution = 1)

plot1 <- DimPlot(merged_seurat_filtered, reduction = "tsne", label = F)
plot2 <- DimPlot(merged_seurat_filtered, reduction = "umap", label = F)

plot1 + plot2

# by grouping examining the batch effect
plot3 <- DimPlot(merged_seurat_filtered, reduction = "tsne", label = F, group.by = 'orig.ident')
plot4 <- DimPlot(merged_seurat_filtered, reduction = "umap", label = F, group.by ='orig.ident')

plot3 + plot4



# Now integrate data

# Integrate data to correct batch effects

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

plot1 + plot2

# by grouping examining the batch effect
plot3 <- DimPlot(seurat.integrated, reduction = "tsne", label = F, group.by = 'orig.ident')
plot4 <- DimPlot(seurat.integrated, reduction = "umap", label = F, group.by ='orig.ident')
plot3+plot4