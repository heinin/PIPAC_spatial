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

data <- read_h5ad("/scratch/hnatri/PIPAC/cell_merged_spatial_filtered_splitsamples_clustered_NC50_NN20_PC20_AnnData.h5ad")
seurat <- CreateSeuratObject(counts = t(as.matrix(data$raw)), 
                             meta.data = data$obs)

head(seurat@meta.data)

#cell_exp <- LayerData(seurat,
#                      assay = "RNA",
#                      layer = "counts")

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

# Adding nuclear counts
merged_spatial <- readRDS("/scratch/hnatri/PIPAC/cell_merged_spatial_filtered_splitsamples.rds")

seurat@assays$nucleus_RNA <- merged_spatial@assays$nucleus_RNA

# Saving as Seurat
saveRDS(seurat, "/scratch/hnatri/PIPAC/cell_merged_spatial_filtered_splitsamples_clustered_NC50_NN20_PC20_Seurat.rds")
