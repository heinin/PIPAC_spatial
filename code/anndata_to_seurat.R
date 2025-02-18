#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/05/2025
# Description: Converting AnnData to Seurat
#==============================================================================#

#==============================================================================
# Libraries
#==============================================================================

library(Seurat)
library(reticulate)
library(anndata)

#==============================================================================
# Convert data
#==============================================================================

data <- read_h5ad("/scratch/hnatri/PIPAC/nonimmune_complete_subset_clustered_NC50_NN20_PC20_AnnData.h5ad")
seurat <- CreateSeuratObject(counts = t(as.matrix(data$X)), 
                             meta.data = data$obs)

head(seurat@meta.data)

# Add PCA reduction
pca <- data$obsm$X_pca
colnames(pca) <- paste0("PCA_", 1:50)
rownames(pca) <- colnames(seurat)
seurat@reductions[["pca"]] <- CreateDimReducObject(embeddings = pca, key = "PCA_", assay = DefaultAssay(seurat))

# Add spatial coordinates
position_xy <- cbind((seurat$x_centroid)*-1,
                     (seurat$y_centroid))
row.names(position_xy) <- row.names(seurat@meta.data)
colnames(position_xy) <- c("SP_1", "SP_2")
seurat[["sp"]] <- CreateDimReducObject(
  embeddings = position_xy, key = "SP_", assay = DefaultAssay(seurat))

#spatial <- data$obsm$SP
#colnames(spatial) <- paste0("SP_", 1:2)
#rownames(spatial) <- colnames(seurat)
#seurat@reductions[["sp"]] <- CreateDimReducObject(embeddings = spatial, key = "SP_", assay = DefaultAssay(seurat))

DimPlot(seurat,
        group.by = "TMA",
        split.by = "TMA",
        ncol = 3,
        #cols = cluster_col,
        reduction = "sp",
        raster = T,
        label = F) +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  NoLegend() +
  RotatedAxis()

# Add UMAP reduction
umap <- data$obsm$X_umap
colnames(umap) <- paste0("UMAP_", 1:2)
rownames(umap) <- colnames(seurat)
seurat@reductions[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = DefaultAssay(seurat))

DimPlot(seurat,
        reduction = "umap",
        group.by = "Sample",
        raster = T) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  NoLegend() +
  ggtitle("")

# Saving as Seurat
saveRDS(seurat, "/scratch/hnatri/PIPAC/nonimmune_complete_subset_clustered_NC50_NN20_PC20_Seurat.rds")
