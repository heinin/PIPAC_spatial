#==============================================================================#
# Author(s) : Heini M Natri, hnatri@tgen.org
# Date: 02/05/2025
# Description: Converting AnnData to Seurat
#==============================================================================#

#==============================================================================
# Libraries
#==============================================================================

library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(reticulate)
library(anndata)
library(ggplot2)

#==============================================================================
# Convert data
#==============================================================================

seurat_data <- readRDS("/scratch/hnatri/PIPAC/stromal_seurat_data.rds")

DimPlot(seurat_data,
        group.by = "leiden_1.5",
        #cols = pipac_cell_immune_clusters_col,
        reduction = "umap",
        raster = T,
        label = T) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  NoLegend()

# Importing each element and building the seurat object
raw_counts <- read.csv("/scratch/hnatri/PIPAC/RSC_latest_EDM_2025-08-06/stromal_raw_counts.csv", header = F)
obs <- read.csv("/scratch/hnatri/PIPAC/RSC_latest_EDM_2025-08-06/stromal_obs.csv")
pcs <- as.matrix(read.csv("/scratch/hnatri/PIPAC/RSC_latest_EDM_2025-08-06/stromal_pcs.csv", header = F))
umap <- as.matrix(read.csv("/scratch/hnatri/PIPAC/RSC_latest_EDM_2025-08-06/stromal_umap.csv", header = F))

rownames(obs) <- obs$cell_id
obs <- obs[,-which(names(obs) %in% c("X"))]

colnames(raw_counts) <- rownames(seurat_data)
rownames(raw_counts) <- obs$cell_id
rownames(pcs) <- colnames(seurat_data)
rownames(umap) <- colnames(seurat_data)

seurat <- CreateSeuratObject(counts = t(as.matrix(raw_counts)), 
                             meta.data = obs)

# Add PCA reduction
colnames(pcs) <- paste0("PCA_", 1:100)
seurat@reductions[["pca"]] <- CreateDimReducObject(embeddings = pcs, key = "PCA_", assay = DefaultAssay(seurat))

# Add UMAP reduction
colnames(umap) <- paste0("UMAP_", 1:2)
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
#seurat@assays$nucleus_RNA <- seurat_data@assays$nucleus_RNA
#
#head(seurat@assays$RNA$counts)
#head(seurat@assays$nucleus_RNA$counts)

seurat_data@reductions$pca <- seurat@reductions$pca
seurat_data@reductions$umap <- seurat@reductions$umap

seurat_data$leiden_0.5 <- seurat$leiden_0.5
seurat_data$leiden_1.0 <- seurat$leiden_1.0
seurat_data$leiden_1.5 <- seurat$leiden_1.5
seurat_data$leiden_2.0 <- seurat$leiden_2.0

# Saving as Seurat
saveRDS(seurat_data, "/tgen_labs/banovich/PIPAC/Seurat/stromal_clustered_NN20_PC20_Seurat.rds")
