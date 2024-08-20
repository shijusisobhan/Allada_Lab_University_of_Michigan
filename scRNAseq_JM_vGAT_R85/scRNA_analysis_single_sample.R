
library(Seurat)

# Create a Seurat object

counts <- Read10X(data.dir = "X:/Sequencing_data/scRNA_seq_JM_8_14-2024/scRNA_11358-JM/10x_analysis_11358-JM/Sample_11358-JM-1/filtered_feature_bc_matrix/")
seurat <- CreateSeuratObject(counts, project="vGAT_ZT0_1")

# Quality control

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Exclude the cell with
 #1. fewer than 1000 or more than 6000 genes
# 2. fewer than 6000 or more than 75000 total UMI
# 3. mt percentage >5%

seurat <- subset(seurat, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & 
                   nCount_RNA > 6000 & nCount_RNA < 75000 &
                   percent.mt < 5)

VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalization
seurat <- NormalizeData(seurat)

# Feature selection for following heterogeneity analysis

seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
top_features <- head(VariableFeatures(seurat), 10)
plot1 <- VariableFeaturePlot(seurat)
LabelPoints(plot = plot1, points = top_features, repel = TRUE)


# Data scaling
seurat <- ScaleData(seurat)
seurat <- ScaleData(seurat, vars.to.regress = c("nFeature_RNA", "percent.mt"))

# PCA

seurat <- RunPCA(seurat, npcs = 50)
ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))
# Plot the cells in the 2D PCA projection
DimPlot(seurat, reduction = "pca")

# Non-linear dimension reduction for visualization
seurat <- RunTSNE(seurat, dims = 1:20)
seurat <- RunUMAP(seurat, dims = 1:20)

plot1 <- TSNEPlot(seurat)
plot2 <- UMAPPlot(seurat)
plot1 + plot2

# Cluster the cells
seurat <- FindNeighbors(seurat, dims = 1:20)
seurat <- FindClusters(seurat, resolution = 1)

plot1 <- DimPlot(seurat, reduction = "tsne", label = TRUE)
plot2 <- DimPlot(seurat, reduction = "umap", label = TRUE)
plot1 + plot2